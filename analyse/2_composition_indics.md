# Indicators at the composition level

Add the columns for the attributes
```sql
ALTER TABLE bd ADD COLUMN cp_fac_lgt REAL DEFAULT 0.0,
               ADD COLUMN cp_fac_lgt_rel REAL DEFAULT 0.0,

               ADD COLUMN cp_d_nxtbuil REAL,
               ADD COLUMN cp_d_builnxtblk REAL,
               ADD COLUMN cp_infzn_area REAL,
               ADD COLUMN cp_r_ftp_infzon REAL,

               ADD COLUMN cp_d_blklim REAL,
               ADD COLUMN cp_blk_area REAL,
               ADD COLUMN cp_r_ftp_blk REAL,
               ADD COLUMN cp_permeab INT DEFAULT 0,
               ADD COLUMN cp_closeness REAL DEFAULT 0.0;
```

## Create _VIEWS_ for the recurring tables
```sql
WITH agglo_buff AS (
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
), ag_ext AS (
  SELECT ST_SetSRID(ST_Extent(geom), 21781) AS geom
  FROM agglo_buff
)
SELECT * INTO ag_ext FROM ag_ext;
, built AS (
  SELECT ST_Union(geom) AS geom
  FROM bd
), non_built AS (
  SELECT ST_Difference(ag_ext.geom, built.geom) AS geom
  FROM ag_ext, built
)
SELECT * INTO toutbatid FROM non_built;
```
And the query to compute all influence zones in R. This takes a lot of time to be computed
```R
#install.packages('rpostgis', 'magrittr')
library(rpostgis)
library(magrittr)

conn <- dbConnect("PostgreSQL", dbname = 'base_data', host = 'localhost', user = 'postgres', password = 'asd')

res <- dbGetQuery(conn, "WITH agglo_buff AS (
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
)
SELECT box2d(geom) FROM agglo_buff;")

bbox <- gsub(")", '', res) %>%
    gsub("BOX(", '', ., fixed = T) %>%
    gsub(',', ' ', .)

bbox <- as.numeric(unlist(strsplit(bbox, ' ')))
xmin <- round(bbox[1], -4)
xmax <- bbox[3]
ymin <- round(bbox[2], -4)
ymax <- bbox[4]

initial_query <- "CREATE TABLE infzone (
  id INT,
  area DOUBLE PRECISION
);"

rec_query <- "WITH rect AS ( -- Rectangle
  SELECT ST_GeomFromText('POLYGON ((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781) AS geom
), extent AS ( -- Prendre un buffer de 550m
  SELECT ST_Buffer(geom, 550) AS geom
  FROM rect
), bats AS ( -- Selectionner les bats dans la zone
  SELECT id, a.geom
  FROM bd a, extent b
  WHERE ST_ContainsProperly(b.geom, a.geom)
), norm_bat AS (
  SELECT *
  FROM bats
  WHERE id NOT IN (SELECT a.id FROM bats a, bats b WHERE ST_Intersects(a.geom, b.geom) AND a.id <> b.id)
), built AS (
  SELECT ST_Union(geom) AS geom
  FROM norm_bat
), non_built AS (
  SELECT ST_Difference(a.geom, b.geom) AS geom
  FROM extent a, built b
), lim_inf_zone AS ( -- bords de la zone d'influence
  SELECT ST_ApproximateMedialAxis(geom) AS geom
  FROM non_built
), inf_zone AS ( -- Zone d'influence
  SELECT (ST_Dump(ST_Polygonize(geom))).geom AS geom
  FROM lim_inf_zone
), bat_ar AS ( -- aire de la zone du batiment
  SELECT id, ST_Area(b.geom) AS area
  FROM bd a, inf_zone b
  WHERE ST_Intersects(a.geom, b.geom)
), bats_def AS ( -- Join the info to the geom
  SELECT *
  FROM bats
  NATURAL INNER JOIN bat_ar  
), bats_rect AS ( -- keep only the building in the zone of interest
  SELECT a.*
  FROM bats_def a, rect b
  WHERE ST_Intersects(a.geom, b.geom)
)
INSERT INTO infzone(id, area)
SELECT id, area
FROM bats_rect;"

sub_in_query <- function(xmin, xmax, ymin, ymax, query) {
    q <- gsub('XMIN', xmin, query) %>%
         gsub('XMAX', xmax, .) %>%
         gsub('YMIN', ymin, .) %>%
         gsub('YMAX', ymax, .)
    return(q)
}

dbSendQuery(conn, initial_query)
for(x in seq(xmin, xmax + 1499, by=1500)) {
    start_x <- x
    end_x <- x + 1500
    for(y in seq(ymin, ymax + 1499, by=1500)) {
        start_y <- y
        end_y <- y + 1500

        q <- sub_in_query(start_x, end_x, start_y, end_y, rec_query)
        dbSendQuery(conn, q)
    }
}

permeab_query <- "WITH pts_blk AS (
  SELECT blk_id, ST_Union(ST_PointOnSurface(geom)) AS geom
  FROM (SELECT id, blk_id, geom FROM bd WHERE blk_id = 9710) z
  GROUP BY blk_id
), vorpol AS (
  SELECT blk_id, ST_CollectionExtract(ST_VoronoiPolygons(geom), 3) AS geom
  FROM pts_blk
), vorlines AS (
  SELECT blk_id, ST_ExteriorRing((ST_Dump(geom)).geom) AS geom
  FROM vorpol
), vorlinesbis AS (
  SELECT blk_id, ST_Union(geom) AS geom
  FROM vorlines
  GROUP BY blk_id
), blk_bord AS (
  SELECT bid, ST_ExteriorRing((ST_Dump(a.geom)).geom) AS geom
  FROM bl a, vorlinesbis b
  WHERE bid = blk_id
), intersec AS (
  SELECT blk_id, (ST_Dump(ST_Intersection(a.geom, b.geom))).geom AS geom
  FROM vorlinesbis a, blk_bord b
  WHERE ST_Intersects(a.geom, b.geom)
  AND a.blk_id = b.bid
), int_count AS (
  SELECT blk_id, COUNT(*) AS tot
  FROM intersec
  GROUP BY blk_id
)
UPDATE bd SET cp_permeab = c.tot / ST_Perimeter(b.geom)
             FROM bd a, bl b, int_count c
             WHERE a.blk_id = b.bid
             AND a.blk_id = c.blk_id
             AND b.bid = c.blk_id;"

ids <- t(dbGetQuery(conn, "SELECT DISTINCT blk_id FROM bd ORDER BY blk_id;"))

for(id in ids) {

}

```

### Compute the street related indicator
```sql
WITH buff_blk AS ( -- Negative buffer on block
  SELECT bid, (ST_Dump(ST_Buffer(geom, -3))).geom AS geom
  FROM bl
), lim_buff AS ( -- Get the boundary of the buffer
  SELECT bid, ST_ExteriorRing(geom) AS geom
  FROM buff_blk
), intersec AS ( -- Get the intersection of the building and buffer's boundary
  SELECT id, ST_Intersection(a.geom, b.geom) AS geom
  FROM bd a, lim_buff b
  WHERE ST_Intersects(a.geom, b.geom)
), build_out_buff AS ( -- Get building's parts which are out of the buffer
  SELECT a.id, (ST_Dump(ST_Difference(a.geom, b.geom))).geom AS geom
  FROM bd a, buff_blk b
  WHERE ST_Intersects(a.geom, b.geom)
), out_length AS ( -- Get the total length of the outer parts
  SELECT id, SUM(ST_Length(ST_ExteriorRing(geom))) AS o_len
  FROM build_out_buff
  GROUP BY id
), int_length AS ( -- Get the length of the intersection (inside building)
  SELECT id, ST_Length(geom) AS i_len
  FROM intersec
), tot_length AS ( -- Length of facade = outer length - intersection with building
  SELECT a.id, o_len - i_len AS fac_len
  FROM out_length a, int_length b
  WHERE a.id = b.id
)
UPDATE bd SET cp_fac_lgt = fac_len
          FROM tot_length
          WHERE bd.id = tot_length.id;

UPDATE bd SET cp_fac_lgt_rel = cp_fac_lgt / b_perim;
```

### Compute the building related indicators
```sql
-- Distance to nearest building
WITH mindist AS (  
  SELECT a.id, c.dmin
  FROM bd a
  CROSS JOIN LATERAL
  (SELECT
    id,
    a.geom <-> b.geom AS dmin
    FROM bd b
    WHERE a.id <> b.id
    ORDER BY a.geom <-> b.geom
    LIMIT 1
  ) AS c
)
UPDATE bd SET cp_d_nxtbuil = dmin
             FROM mindist
             WHERE bd.id = mindist.id;

-- Distance to nearest building in another block
-- Takes time >505000 ms
WITH mdist_nxt_blk AS (  
 SELECT a.id, c.dmin
 FROM bd a
 CROSS JOIN LATERAL
 (SELECT
   id, blk_id,
   a.geom <-> b.geom AS dmin
   FROM bd b
   WHERE a.blk_id <> b.blk_id
   ORDER BY a.geom <-> b.geom
   LIMIT 1
 ) AS c
)
UPDATE bd SET cp_d_builnxtblk = dmin
          FROM mdist_nxt_blk
          WHERE bd.id = mdist_nxt_blk.id;

-- Influence zone area
WITH uniq_inf_zone AS (
  SELECT id, MAX(area) AS area
  FROM infzone
  GROUP BY id
), inf_zone_max AS (
  SELECT id, CASE WHEN area > 282743.3 THEN 282743.3 ELSE area END AS area
  FROM uniq_inf_zone
)
UPDATE bd SET cp_infzn_area = area
          FROM inf_zone_max
          WHERE bd.id = inf_zone_max.id;

-- Ratio building footprint / influence zone
UPDATE bd SET cp_r_ftp_infzon = b_area / cp_infzn_area;
```

### Compute the block related indicators

```sql
-- Distance to block limit
-- Takes time
WITH blk_bord AS (
  SELECT bid, ST_ExteriorRing((ST_Dump(geom)).geom) AS geom
  FROM bl
), mdist_bord AS (  
 SELECT a.id, c.dmin
 FROM bd a
 CROSS JOIN LATERAL
 (SELECT
   id, blk_id,
   a.geom <-> b.geom AS dmin
   FROM blk_bord b
   WHERE a.blk_id = b.bid
   ORDER BY a.geom <-> b.geom
   LIMIT 1
 ) AS c
)
UPDATE bd SET cp_d_blklim = dmin
          FROM mdist_bord
          WHERE bd.id = mdist_bord.id;

-- Area of block
UPDATE bd SET cp_blk_area = ST_Area(bl.geom)
          FROM bl
          WHERE bd.blk_id = bl.bid;

-- Ratio building footprint to block area
UPDATE bd SET cp_r_ftp_blk = b_area / cp_blk_area;

-- Block permeability
WITH pts_blk AS (
  SELECT blk_id, ST_Union(ST_PointOnSurface(geom)) AS geom
  FROM (SELECT id, blk_id, geom FROM bd WHERE id BETWEEN 1 AND 100) z
  GROUP BY blk_id
), vorpol AS (
  SELECT blk_id, ST_CollectionExtract(ST_VoronoiPolygons(geom), 3) AS geom
  FROM pts_blk
), vorlines AS (
  SELECT blk_id, ST_ExteriorRing((ST_Dump(geom)).geom) AS geom
  FROM vorpol
), vorlinesbis AS (
  SELECT blk_id, ST_Union(geom) AS geom
  FROM vorlines
  GROUP BY blk_id
), blk_bord AS (
  SELECT bid, ST_ExteriorRing((ST_Dump(geom)).geom) AS geom
  FROM bl
  WHERE bid IN (SELECT blk_id FROM bd)
), intersec AS (
  SELECT blk_id, (ST_Dump(ST_Intersection(a.geom, b.geom))).geom AS geom
  FROM vorlinesbis a, blk_bord b
  WHERE ST_Intersects(a.geom, b.geom)
  AND a.blk_id = b.bid
), int_count AS (
  SELECT blk_id, COUNT(*) AS tot
  FROM intersec
  GROUP BY blk_id
)
UPDATE bd SET cp_permeab = c.tot / ST_Perimeter(b.geom)
             FROM bd a, bl b, int_count c
             WHERE a.blk_id = b.bid
             AND a.blk_id = c.blk_id
             AND b.bid = c.blk_id;

-- Block closeness
WITH bloc_perim AS (
 SELECT bid, ST_Perimeter(geom) AS blperim
 FROM bl
 WHERE bid IN (SELECT blk_id FROM bd WHERE cp_fac_lgt <> 0)
), csn AS (
  SELECT a.id, a.cp_fac_lgt / b.blperim AS csns
  FROM bd a, bloc_perim b
  WHERE a.blk_id = b.bid
)
UPDATE bd SET cp_closeness = csns
          FROM csn
          WHERE bd.id = csn.id;
```

### Miscellaneous
SQL query for the R script
```sql
-- TEST 3
WITH rect AS ( -- Rectangle
  SELECT ST_GeomFromText('POLYGON ((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781) AS geom
), extent AS ( -- Prendre un buffer de 550m
  SELECT ST_Buffer(geom, 550) AS geom
  FROM rect
), bats AS ( -- Selectionner les bats dans la zone
  SELECT id, a.geom
  FROM bd a, extent b
  WHERE ST_ContainsProperly(b.geom, a.geom)
), norm_bat AS (
  SELECT *
  FROM bats
  WHERE id NOT IN (SELECT a.id FROM bats a, bats b WHERE ST_Intersects(a.geom, b.geom) AND a.id <> b.id)
), built AS (
  SELECT ST_Union(geom) AS geom
  FROM norm_bat
), non_built AS (
  SELECT ST_Difference(a.geom, b.geom) AS geom
  FROM extent a, built b
), lim_inf_zone AS ( -- bords de la zone d'influence
  SELECT ST_ApproximateMedialAxis(geom) AS geom
  FROM non_built
), inf_zone AS ( -- Zone d'influence
  SELECT (ST_Dump(ST_Polygonize(geom))).geom AS geom
  FROM lim_inf_zone
), bat_ar AS ( -- aire de la zone du batiment
  SELECT id, ST_Area(b.geom) AS area
  FROM bd a, inf_zone b
  WHERE ST_Intersects(a.geom, b.geom)
), bats_def AS ( -- Join the info to the geom
  SELECT *
  FROM bats
  NATURAL INNER JOIN bat_ar  
), bats_rect AS ( -- keep only the building in the zone of interest
  SELECT a.*
  FROM bats_def a, rect b
  WHERE ST_Intersects(a.geom, b.geom)
)
INSERT INTO infzone(id, area)
SELECT id, area
FROM bats_rect;
```
