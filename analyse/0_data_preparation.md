# Data preprocessing

First of all, prepare the DB:
```sql
CREATE DATABASE base_data;
-- Don't forget to connect to the DB \c base_data
CREATE EXTENSION postgis_sfcgal CASCADE;
CREATE EXTENSION pgrouting;
```


## Workaroud the SWISSBUILDINGS3D shapefile
A python script (coded by C. Kaiser) is used to extract the simple geometries from the features and output a sql file to insert them in the DB:
```python
import ogr
import json

path = "D:\Memoire\data\swiss_bldg"

def SimpleGeometryFromFeature(feat):
    geom = feat.geometry()
    max_area = -1.0
    max_area_geom = None
    geom_cnt = geom.GetGeometryCount()
    for i in range(geom_cnt):
        gi = geom.GetGeometryRef(i)
        gi_area = gi.Area()
        if gi_area > max_area:
            max_area = gi_area
            max_area_geom = gi
    return max_area_geom


def AttributesFromFeature(feat, attrs):
    feat_json = json.loads(feat.ExportToJson())
    props = feat_json['properties']
    return [str(props[attr]) for attr in attrs]


shp = "Buildings3D.shp"
driver = ogr.GetDriverByName('ESRI Shapefile')
ds = driver.Open(path, 0)
lyr = ds.GetLayer()
fc = lyr.GetFeatureCount()


# Outfile in SQL format with the WKT geometries and all attributes.
fout = open('buildings3d.sql', 'w')

# Create the table. We assume all attributes are strings.
# Conversion is possible in PostgreSQL with casting (e.g. COLNAME::integer)
feat = lyr.GetFeature(0)
feat_json = json.loads(feat.ExportToJson())
h = feat_json['properties'].keys()

fout.write("CREATE TABLE bldg3d ( fid serial primary key, \n")
fout.write("".join(["%s varchar(250),\n" % col for col in h]))
fout.write("geom GEOMETRY );\n\n")

cols = list(h) + ['geom']

# And now write all the features
for feat in lyr:
    geom = SimpleGeometryFromFeature(feat)
    attributes = AttributesFromFeature(feat, h)
    attributes.append(geom.ExportToWkt())
    fout.write("INSERT INTO buildings3d (%s) VALUES ('%s');\n" % (",".join(cols), "','".join(attributes)))

fout.close()
```

Then import the file in the DB :
`psql -U postgres -d base_data -f buildings3d.sql`

In the DB, keep only the buildings inside a smaller area instead of the whole of Switzerland:
```sql
-- Select only buildings in studied part of Switzerland
SELECT * INTO bldg3d FROM buildings3d
WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((494000 135000, 494000 181000, 572000 181000, 572000 135000, 494000 135000))'));
-- Drop the previous table with all buildings
DROP TABLE buildings3d;
-- Set SRID to 21781
  -- Transform everything into polygons
UPDATE bldg3d SET geom = (ST_Dump(geom)).geom;
  -- Drop the Z index
ALTER TABLE bldg3d
  ALTER COLUMN geom TYPE geometry(Polygon)
    USING ST_Force2d(geom);
  -- Assign 21781 as SRID
SELECT UpdateGeometrySRID('bldg3d', 'geom', '21781');
  -- Add index to the table
CREATE INDEX bldg3d_gix ON bldg3d USING GIST (geom);
VACUUM ANALYZE bldg3d;
```


## Reducing the area of interest
To ensure a lighter database and some faster process, the base datasets are cut according to a bounding box in order to import fewer elements in the DB.
```
:: Cut the original files
ogr2ogr -f "ESRI Shapefile" tlm_roads.shp TLM_STRASSE.shp -clipsrc 494000 135000 572000 181000

ogr2ogr -f "ESRI Shapefile" tlm_bldgs.shp TLM_GEBAEUDE_FOOTPRINT.shp -clipsrc 494000 135000 572000 181000
```

## Inserting into DB

Then prepare and import the shapefiles:
```
:: Prepare the sql file
shp2pgsql -I -s 21781 ./tlm_bldgs.shp bldgs > importBldgs.sql

shp2pgsql -I -s 21781 ./tlm_roads.shp roads > importRoads.sql

:: Import it
psql -U postgres -d base_data -f importBldgs.sql
psql -U postgres -d base_data -f importRoads.sql
```

## Buildings preparation
The TLM doen't include information about buildings' heights, which the SWISSBUILDINGS3D 1.0 has. Therefore we assign the height of the latter to the former where the buildings intersect:
```sql
-- Add a column to be populated
ALTER TABLE bldgs ADD COLUMN b_height REAL;
-- Check the heights of intersecting buildings
WITH int_hei AS (
  SELECT gid, b.height::real AS hei
  FROM bldgs a, bldg3d b
  WHERE ST_Intersects(a.geom, b.geom)
)
UPDATE bldgs SET b_height = hei
             FROM int_hei
             WHERE bldgs.gid = int_hei.gid;
-- UPDATE 120300 lignes sur 129204

DELETE FROM bldgs
  WHERE b_height IS NULL
  OR b_height <= 0;
```

## Roads preparation
For further analyses, we need to assign a width to each road. We take the type (_art_) and match it to a width. A conservative approach is taken in order to avoid roads crossing through buildings:

```sql
-- Add width info to roads
ALTER TABLE roads ADD COLUMN tru_width REAL;
UPDATE roads SET tru_width = CASE
                                WHEN objektart = 'Autobahn' OR objektart = 'Autostrasse' THEN 10.22
                                WHEN objektart = '10m Strasse' THEN 10.21/2
                                WHEN objektart = '6m Strasse' THEN 6.21/2
                                WHEN objektart = '4m Strasse' THEN 4.21/2
                                WHEN objektart = '3m Strasse' THEN 2.81/2
                                WHEN objektart = '2m Weg' OR objektart = '2m Wegfragment' THEN 1.81/2
                                WHEN objektart = '1m Weg' OR objektart = '1m Wegfragment' THEN 0.91/2
                                ELSE 2.81/2
                              END;
```

## Define the agglomerations

### Statistical - OFS
The FSO retains several categories for the communes comprised in an agglomeration :
- Center city
- Main center
- Suburbs

Two statistical agglomerations are retained, the first is limited to the main center around the city (statA), the second englobes the suburbs (statB) and is therefore the largest definition.

Import the file with the communes comprised in the agglomeration and their types:
```sql
CREATE TABLE comms (
  id int, name varchar(100), ct varchar(2), cat varchar(20), pop12 int
);

COPY comms FROM 'D:\Memoire\data\agglos\OFS\agglo_Lausanne.csv' DELIMITER ';' CSV HEADER;
```
Make the same with the shapefile of the communes:
```
shp2pgsql -I -s 21781 ./swissBOUNDARIES3D_1_1_TLM_HOHEITSGEBIET.shp comms_geo > importcomms.sql
psql -U postgres -d base_data -f importcomms.sql
```

##### - StatA
```sql
-- Extract the area of interest
WITH comms_cent AS (-- Select comms in center
  SELECT *
  FROM comms
  WHERE cat IN ('Ville-centre', 'centre_principal')
), com_cent_geo AS (-- Select their geom
  SELECT gid, a.name, geom
  FROM comms_geo a, comms_cent b
  WHERE a.bfs_nummer = b.id
), com_uni AS (-- Select strict area of study
  SELECT ST_Union(geom) AS geom
  FROM com_cent_geo
)
SELECT * INTO stata FROM com_uni;
```

##### - StatB
```sql
WITH comm_agglo AS (
  SELECT gid, a.name, geom
  FROM comms_geo a, comms b
  WHERE a.bfs_nummer = b.id
), agglo_uni AS (
  SELECT ST_Union(geom) AS geom
  FROM comm_agglo
)
SELECT * INTO statb FROM agglo_uni;
```

### Morphological - Fractals
To determine the morphological agglomeration, _Morpholim 1.5_ is used. First we use the tlm_bldg shapefile and let the program run. Then we select the polynomial best describing the curve, in this case 6th. This approximation is used to determine the points of maximum curvature, corresponding to the threshold defining the urban enveloppe.

The urban enveloppes are computed and we then import the shapefile in the DB
```
shp2pgsql -I -s 21781 ./tlm_bldgs-env-89.10986982180968.shp urbenv > importMorphFrac.sql

psql -U postgres -d base_data -f importMorphFrac.sql
```
We then select the main cluster as the morphological agglomeration:
```sql
WITH morph_agglo AS (
  SELECT geom
  FROM urbenv
  ORDER BY ST_Area(geom) DESC
  LIMIT 1
)
SELECT * INTO morfrac FROM morph_agglo;
```

### Morphological - City clustering algorithm (CCA)
To determine the morphological agglomeration using the CCA, the _R_ package `osc` is used.
After a series of test, the threshold of at least 20 inhabitants and/or employments has been retained. This threshold was used on the whole of Switzerland to find the type of neighbourhood which gives the best results.

Two _tolerances_ (the maximal distance separating hectares accepted to consider them adjacent) have been retained :
- king neighbourhood, 142 meters : allows only direct neighbours, included the ones in diagonal (touching by the corners)
- lenient, 223 meters : allows for a gap of maximum one hectare diagonally

```R
# Prepare the environment
needed_packages <- c('osc', 'tidyverse', 'raster', 'rgeos', 'rgdal', 'maptools', 'rpostgis')
# install.packages(needed_packages)
lapply(needed_packages, library, character.only = T)

# Import and prepare the data
statpop <- read.csv("./STATPOP2016G.csv", header = T)
statemp <- read.csv("./STATENT2015_N08_V170824G.csv", header = T)

statpop <- statpop[, c(1:3, 6)]
statemp <- statemp[, c(3:5, 10)]
statemp <- statemp[, c("RELI", "X_KOORD", "Y_KOORD", "B1508EMPT")]

statemp[, 2:3] <- statemp[, 2:3] + 50
statpop[, 2:3] <- statpop[, 2:3] + 50

# Return all hectares with info and sum if inhabitants + employments
allHects <- function(r) {
    if(is.na(r[2])) {
        r[2] <- r[5]
        r[3] <- r[6]
        r[4] <- r[7]
        return(r[1:4])
    } else if(!(is.na(r[7]))) {
        r[4] <- r[4] + r[7]
        return(r[1:4])
    } else {
        return(r[1:4])
    }
}

# Compute it (optionnaly write the result to a file)
hects <- merge(statpop, statemp, by='RELI', all=T)
hects <- apply(hects, 1, function(x) allHects(x))
hects <- data.frame(t(hects))
colnames(hects)[4] <- 'TOT'
# write.csv(hects, file = 'occupiedHects.csv', row.names = F)

# Transform the data into a raster
swissCRS <- CRS("+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
hects_r <- hects[,2:4] %>%
  rasterFromXYZ()
proj4string(hects_r) <- swissCRS
#plot(hects_r)

# Compute the CCA according to the defined tolerance, write the unioned result to the DB
# Outputs a plot of the results at the swiss scale
applyCCA <- function(tolerance, file_name, placeholder) {
    # Perform the CCA
    cca_rslt <- cca(hects_r, cell.class = c(20:max(hects$TOT)), s = tolerance, unit = 100)

    # Separate the results to clarify
    clust_coord <- cca_rslt[[1]]
    clust_size <- cca_rslt[[2]]

    # Refine the analysis to keep only clusters > 3km^2  (300 hects)
    big_enough <- integer()
    for(i in 1:length(clust_size)) {
        if(clust_size[i] > 300) {
            big_enough <- c(big_enough, i)
        }
    }
    placeholder <- clust_coord[clust_coord[, 3] %in% big_enough, ]
    rownames(placeholder) <- NULL

    # Plot the resulting clusters
    fig <- ggplot(placeholder, aes(x = long, y = lat)) +
              geom_point(shape = 15, size = 0.05, color = placeholder$cluster_id)
    leg <- labs(title = "City clusters of Switzerland", subtitle = paste('>20 inhabitants and/or employments, tolerance', tolerance, "m", sep = ' '), x='', y='', caption = 'Sources : OFS, site web Statistique suisse, 2017. Package osc.')
    # Transform the result as a SpatialPointsDF
    coordinates(placeholder) <- 1:2
    proj4string(placeholder) <- swissCRS

    # Make a square buffer of 50m to cover the hectare
    hects_poly <- gBuffer(placeholder, byid = T, width = 50, capStyle = "SQUARE")

    # Get the union of the buffer to keep only one polygon per cluster
    clusters <- gUnaryUnion(hects_poly, id = hects_poly$cluster_id)
    IDs <- data.frame(ID = sapply(slot(clusters, "polygons"), function(x) slot(x, "ID")))
    rownames(IDs)  <- IDs$ID
    clusters <- SpatialPolygonsDataFrame(clusters, IDs)
    proj4string(clusters) <- swissCRS

    # Export the result in Postgis
    conn <- dbConnect("PostgreSQL", dbname = 'base_data', host = 'localhost', user = 'postgres', password = 'asd')
    pgInsert(conn, file_name, clusters, new.id = "gid")
    # Adjust SRID in DB
    file_name <- paste("'", file_name, "'", sep = '')
    q <- paste("SELECT UpdateGeometrySRID(", file_name, ",'geom', 21781)", sep = '')
    dbSendQuery(conn, q)

    # Return plot
    return(fig + leg)                            
}

# CCA with a king neighbourhood
applyCCA(142, 'cca_king', 'cca')

# CCA permitting a gap of more than one hectare
applyCCA(224, 'cca_lenient', 'cca2')
```

Since all the agglomerations of Switzerland are added to the DB, it must be refined to the main agglomeration around Lausanne for both types.

```sql
-- For the CCA_King there are two sufficiently big cluster, so the biggest is conserved
WITH agglos_ls AS (
  SELECT a.*
  FROM cca_k a, statb b
  WHERE ST_Intersects(a.geom, b.geom)
), main_agglo AS (
  SELECT *, ST_Area(geom) AS area
  FROM agglos_ls
  ORDER BY area DESC LIMIT 1
)
SELECT *
INTO cca_king
FROM main_agglo;

DROP TABLE cca_k;

-- For the CCA_Lenient, there is only one cluster, so it is taken as is.
SELECT a.*
INTO cca_lenient
FROM cca_l a, statb b
WHERE ST_Intersects(a.geom, b.geom);

DROP TABLE cca_l;
```

## Roads and buildings selections

The various indicators will be computed for all buildings and then the classification will be made on the buildings located in each agglomeration. So in order to avoid an artificial border effect, we buffer the largest agglomeration retained, in this case _statb_.

```sql
WITH agglo_buff AS ( -- Have a margin of 2.5km
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
) -- Keep only the buildings in this area
SELECT a.gid, a.geom, a.b_height AS b_height, a.objektart
INTO bd
FROM bldgs a, agglo_buff b
WHERE ST_Intersects(a.geom, b.geom);

-- Adjust the IDs and create a spatial index
ALTER TABLE bd DROP COLUMN gid;
ALTER TABLE bd ADD COLUMN id SERIAL PRIMARY KEY;
CREATE INDEX bd_gix ON bd USING GIST (geom);
VACUUM ANALYZE bd;

-- Take the same margin for the roads
WITH agglo_buff AS (
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
)
SELECT a.gid, a.geom, a.tru_width, a.kunstbaute
INTO rd
FROM roads a, agglo_buff b
WHERE ST_Intersects(a.geom, b.geom)
  AND a.kunstbaute IN ('Keine', 'Bruecke', 'Treppe');
```

In order to avoid building next to the lake being considered as bordered only on one side, the lake contour is added in the roads with a tru_width of 0.5m
The contour is issued from the TLM_STEHENDES_GEWAESSER shapefile.

```
shp2pgsql -I -s 21781 ./TLM_STEHENDES_GEWAESSER.shp lakes > importLakes.sql
psql -U postgres -d base_data -f importLakes.sql
```
 Then keep only the contours
```sql

WITH agglo_buff AS (
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
), leman AS ( -- Keep only parts of the lake near the agglo
  SELECT a.geom, a.objektart, a.name
  FROM lakes a, agglo_buff b
  WHERE ST_Intersects(a.geom, b.geom)
  AND a.name = 'Le Léman'
) -- Add them in the roads table
INSERT INTO rd(geom, tru_width, kunstbaute)
SELECT geom, 0.5, 'Keine'
FROM leman;

-- Adjust the IDs and create s spatial index
ALTER TABLE rd DROP COLUMN gid;
ALTER TABLE rd ADD COLUMN id SERIAL PRIMARY KEY;
CREATE INDEX rd_gix ON rd USING GIST (geom);
VACUUM ANALYZE rd;

-- Drop the lakes table
DROP TABLE lakes;
```

## Roads topology
```sql
CREATE EXTENSION postgis_topology;
SET search_path = topology, public;

SELECT topology.CreateTopology('rd_topo', 21781);
SELECT topology.AddTopoGeometryColumn('rd_topo', 'public', 'rd', 'topo_geom', 'LINESTRING');
-- Takes a long time (2143000ms)
UPDATE rd SET topo_geom = topology.toTopoGeom(ST_Force2d(geom), 'rd_topo', 1, 1.0);

-- Use the edges as roads
SELECT edge_id AS id, geom INTO rd_f FROM rd_topo.edge;
-- Add a column to add the width
ALTER TABLE rd_f ADD COLUMN tru_width REAL DEFAULT 2.81/2;
-- Get the width of roads which intersect the edge
WITH edges_int AS (
  SELECT a.edge_id, b.tru_width AS width, ST_Length(ST_Intersection(a.geom, ST_Force2d(b.geom))) AS int_l, a.geom
  FROM rd_topo.edge a, rd b
  WHERE ST_Intersects(a.geom, ST_Force2d(b.geom))
), edges_pos AS ( -- Keep only the main intersection
  SELECT *
  FROM edges_int
  WHERE int_l > 0
)
UPDATE rd_f SET tru_width = width
            FROM edges_pos
            WHERE rd_f.id = edges_pos.edge_id;

```

## Blocks preparation
First the blocks, defined as the space encomprised between roads are added to a table. The road buffer (with their width as defined in the data preparation) are computed and added too.

```sql
WITH blk AS (
  SELECT (ST_Dump(ST_Polygonize(geom))).geom AS geom
  FROM rd_f
), blk_w_bd AS (
  SELECT DISTINCT a.geom
  FROM blk a, bd b
  WHERE ST_Intersects(a.geom, b.geom)
)
SELECT * INTO blk FROM blk_w_bd;
CREATE INDEX blk_gix ON blk USING GIST (geom);
VACUUM ANALYZE blk;

WITH road_buf AS ( -- All the place taken by the roads
  SELECT ST_Buffer(geom, tru_width, 'endcap=square join=round') AS geom
  FROM rd_f
)
SELECT * INTO road_buf FROM road_buf;
CREATE INDEX road_buff_gix ON road_buf USING GIST (geom);
VACUUM ANALYZE road_buf;
```

### SQL IN R
We need to union the results of the road buffers but this would mean a very complex object and some computing time. Therefore we use an R script to compute the result for smaller grids of the area. This method just leave an elongated block empty, thus we add it after, based on its centroid.

```R
### This script compute the blocks minus buffers
### for the area of interest, one grid at a time

## Environment
#install.packages('rpostgis','magrittr')
library(rpostgis)
library(magrittr)

# Create a connection
conn <- dbConnect("PostgreSQL", dbname = 'base_data', host = 'localhost', user = 'postgres', password = 'asd')

# Get the bbox of the area of interest
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

# Query for the first block (creates table)
initial_query <- "SELECT geom
INTO full_blk
FROM blk
WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((510023 144456, 510023 150000, 520000 150000, 520000 144456, 510023 144456))', 21781));"

# Road query
road_query <- "WITH roads_buf AS (
  SELECT geom
  FROM road_buf
  WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781))
), union_buf AS (
  SELECT ST_Union(geom) AS geom
  FROM roads_buf
)"
# Block query (part of the same query)
block_query <- " block AS (
  SELECT geom
  FROM blk
  WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781))
), diff AS (
  SELECT ST_Difference(a.geom, b.geom) AS geom
  FROM block a, union_buf b
)
INSERT INTO full_blk(geom)
SELECT geom
FROM diff;"

# This replace the coordinates in the queries
sub_in_query <- function(xmin, xmax, ymin, ymax, road_q, block_q) {
    if(missing(block_q)) {
        q <- gsub('XMIN', xmin, road_q) %>%
                     gsub('XMAX', xmax, .) %>%
                     gsub('YMIN', ymin, .) %>%
                     gsub('YMAX', ymax, .)
    } else {
        rq <- gsub('XMIN', xmin-2500, road_q) %>%
                     gsub('XMAX', xmax+2500, .) %>%
                     gsub('YMIN', ymin-2500, .) %>%
                     gsub('YMAX', ymax+2500, .)

        bq <- gsub('XMIN', xmin, block_q) %>%
                     gsub('XMAX', xmax, .) %>%
                     gsub('YMIN', ymin, .) %>%
                     gsub('YMAX', ymax, .)
        q <- paste(rq, bq, sep = ',')
    }
    return(q)
}

# Compute the blocks for each square and insert
# into DB
for(x in seq(xmin, xmax + 7499, by=7500)) {
    start_x <- x #- 1000
    end_x <- x + 7500
    for(y in seq(ymin, ymax + 4999, by=5000)) {
        start_y <- y #- 1000
        end_y <- y + 5000
        if(x == xmin && y == ymin) {
            q <- sub_in_query(xmin, end_x, ymin, end_y, initial_query)
            dbSendQuery(conn, q)
        } else {
            q <- sub_in_query(start_x, end_x, start_y, end_y, road_query, block_query)
            dbSendQuery(conn, q)
        }
    }
}
```

Now the blocks are in the table, but there can be some blocks at the border of grid computed twice and some being differentiated only on part of their perimeter. We remove them and add the final block.
```sql
-- Avoid the multiple geometries
SELECT DISTINCT (ST_Dump(geom)).geom AS geom
INTO blocks
FROM full_blk
WHERE st_geometrytype(geom) <> 'ST_GeometryCollection';

-- Select only the inner blocks
SELECT a.geom
INTO bl
FROM blocks a, blk b
WHERE ST_Intersects(a.geom, b.geom)
AND ST_ContainsProperly(b.geom, a.geom);

-- Add the missing block
-- Trade-off of automation for speed
INSERT INTO bl(geom)
SELECT geom
FROM blocks
WHERE ST_GeomFromText('POINT(551506.012223242 165006.729552064)', 21781) = ST_Centroid(blocks.geom);


-- Create a spatial index
ALTER TABLE bl ADD COLUMN bid SERIAL;
CREATE INDEX bl_gix ON bl USING GIST (geom);
VACUUM ANALYZE bl;
```

We must then assign a block to each building. For building spanning across multiple blocks, we take the block with the biggest share of the building.
```sql
ALTER TABLE bd ADD COLUMN blk_id INTEGER;

WITH int_area AS (
  SELECT id, bid, ST_Area(ST_Intersection(a.geom, b.geom)) AS int_a
  FROM bd a, bl b
  WHERE ST_Intersects(a.geom, b.geom)
), max_int AS (
  SELECT a.*
  FROM int_area a, (SELECT id, MAX(int_a) AS max FROM int_area GROUP BY id) b
  WHERE a.id = b.id
    AND a.int_a = b.max
)
UPDATE bd SET blk_id = bid
          FROM max_int
          WHERE bd.id = max_int.id;
```


## OSM
First, we need to have the bounding box of the _statB_ buffer to download the data.

```sql
WITH agglo_buff AS (
  SELECT ST_Buffer(geom, 2500) AS geom
  FROM statb
)
SELECT box2d(ST_Transform(ST_GeometryFromText(ST_AsText(box2d(geom)), 21781), 4326)) FROM agglo_buff;
-- BOX(6.26087181223305 46.4453855400558,6.93584863091059 46.7633604101941)
```

Once the data in the bounding box is downloaded, we can import the .osm file in the DB
```
osm2pgrouting -f statb.osm -c D:\PostgreSQL\9.6\bin\mapconfig.xml -U postgres -W asd -d base_data
```

### Miscellaneous

Queries used in the R script
```sql
-- Initial query
SELECT geom
INTO full_blk
FROM blk
WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((510023 144456, 510023 150000, 520000 150000, 520000 144456, 510023 144456))', 21781));

-- Road, block and difference query
WITH roads_buf AS (
  SELECT geom
  FROM road_buf
  WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781))
), union_buf AS (
  SELECT ST_Union(geom) AS geom
  FROM roads_buf
), block AS (
  SELECT geom
  FROM blk
  WHERE ST_Intersects(geom, ST_GeomFromText('POLYGON((XMIN YMIN, XMIN YMAX, XMAX YMAX, XMAX YMIN, XMIN YMIN))', 21781))
), diff AS (
  SELECT ST_Difference(a.geom, b.geom) AS geom
  FROM block a, union_buf b
)
INSERT INTO full_blk(geom)
SELECT geom
FROM diff;
```
