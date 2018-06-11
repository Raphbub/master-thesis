# Indicators at the building level

Add the columns for the attributes
```sql
ALTER TABLE bd ADD COLUMN b_perim REAL,
               ADD COLUMN b_area REAL,
               ADD COLUMN b_r_vol_fac REAL,
               ADD COLUMN b_MaxEdge REAL,
               ADD COLUMN b_MinEdge REAL,
               ADD COLUMN b_stories INT,
               ADD COLUMN b_floorsqm REAL,

               ADD COLUMN c_Miller REAL,
               ADD COLUMN c_Schumm REAL,
               ADD COLUMN c_Haggett REAL,
               ADD COLUMN c_LeeSallee REAL,
               ADD COLUMN c_Ehrenb REAL,

               ADD COLUMN bb_perim REAL,
               ADD COLUMN bb_area REAL,
               ADD COLUMN bb_length REAL,
               ADD COLUMN bb_width REAL,
               ADD COLUMN bb_r_lw REAL,
               ADD COLUMN bb_r_area REAL,
               ADD COLUMN bb_r_perim REAL,

               ADD COLUMN cc_rad REAL,
               ADD COLUMN cc_exch REAL,
               ADD COLUMN cc_detour REAL,

               ADD COLUMN ch_area REAL,
               ADD COLUMN ch_perim REAL,
               ADD COLUMN ch_r_area REAL,
               ADD COLUMN ch_r_perim REAL,

               ADD COLUMN s_deadend REAL,
               ADD COLUMN sc_lines REAL,
               ADD COLUMN sc_length REAL,
               ADD COLUMN sc_orient INT DEFAULT 0,
               ADD COLUMN sc_l_sn REAL DEFAULT 0.0,
               ADD COLUMN sc_l_ew REAL DEFAULT 0.0,
               ADD COLUMN sc_l_nesw REAL DEFAULT 0.0,
               ADD COLUMN sc_l_senw REAL DEFAULT 0.0,
               ADD COLUMN sc_m_orient VARCHAR(4),

               ADD COLUMN m_corndis REAL,
               ADD COLUMN m_court INT DEFAULT 0,
               ADD COLUMN m_court_area REAL DEFAULT 0.0,
               ADD COLUMN m_court_rel_a REAL DEFAULT 0.0,

               ADD COLUMN dm_inscr_c REAL;
```

## Create _VIEWS_ for the recurring tables and precalculate other needed attributes
```sql
-- View circumscribed circle, with radius extracted
CREATE OR REPLACE VIEW ccirc AS
  SELECT id, (ST_MinimumBoundingRadius(geom)).radius AS rad
  FROM bd;

-- View bounding box
CREATE OR REPLACE VIEW bbox AS
  SELECT id, Box2d(geom) AS bb
  FROM bd;

-- View convex hull
CREATE OR REPLACE VIEW convhull AS
  SELECT id, ST_ConvexHull(geom) AS chull
  FROM bd;

-- View skeleton & centerline
CREATE OR REPLACE VIEW skeleton AS
  SELECT id, ST_StraightSkeleton(geom) AS skel,
         ST_ApproximateMedialAxis(geom) AS ctl
  FROM bd;

-- View for the orientation of the centerlines' segments
CREATE OR REPLACE VIEW sc_orien AS (
  WITH clpts AS (
    SELECT id, (ST_DumpPoints(ctl)).geom AS pts,
           ((ST_DumpPoints(ctl)).path)[1] AS place
    FROM skeleton
  ), aziline AS (
    SELECT a.id, a.place, DEGREES(ST_Azimuth(a.pts, b.pts)) AS orientation, ST_Distance(a.pts, b.pts) AS long
    FROM clpts a, clpts b
    WHERE a.id = b.id
      AND a.place = b.place
      AND NOT a.pts = b.pts
      AND DEGREES(ST_Azimuth(a.pts, b.pts)) BETWEEN 67.5 AND 247.5
  ), orlc AS (
    SELECT *,
      CASE WHEN ROUND(CAST(orientation AS numeric), 2) BETWEEN 67.5 AND 112.5 THEN 'EW'
           WHEN ROUND(CAST(orientation AS numeric), 2) BETWEEN 112.501 AND 157.5 THEN 'SENW'
           WHEN ROUND(CAST(orientation AS numeric), 2) BETWEEN 157.501 AND 202.5 THEN 'SN'
           WHEN ROUND(CAST(orientation AS numeric), 2) BETWEEN 202.501 AND 247.5 THEN 'SWNE'
           ELSE 'ERROR'
      END AS orlc
    FROM aziline
  )
  SELECT * FROM orlc
);

-- Sum of centerlines by orientation
CREATE OR REPLACE VIEW sum_orien AS (
  WITH sum_p_or AS (-- Get the total length by orientation
    SELECT id, orlc, SUM(long) AS ltot
    FROM sc_orien
    GROUP BY id, orlc
  )
  SELECT * FROM sum_p_or
);

-- Computing an approximation of the inscribed circle diameter
WITH ctline AS ( -- Dump centerline
  SELECT id, (ST_Dump(ctl)).geom AS lines,
         (ST_Dump(ctl)).path[1] AS place
  FROM skeleton
), dists AS ( -- Compute distance from line to exterior of polygon
  SELECT bd.id AS id, place, lines,
         ST_Distance(lines, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS dist
  FROM ctline, bd
  WHERE ctline.id = bd.id
), max_dists AS ( -- Find the max distance
  SELECT id, place, lines, dist, MAX(dist) OVER (PARTITION BY id) AS max_dist
  FROM dists
), furthest_lines AS ( -- Find the furthest line
  SELECT *
  FROM max_dists
  WHERE dist = max_dist
), fl_pts AS ( -- Divide the line in 7 points
  SELECT id, ST_StartPoint(lines) AS pta,
         ST_LineInterpolatePoint(lines, 0.2) AS ptb,
         ST_LineInterpolatePoint(lines, 0.4) AS ptc,
         ST_LineInterpolatePoint(lines, 0.5) AS ptd,
         ST_LineInterpolatePoint(lines, 0.6) AS pte,
         ST_LineInterpolatePoint(lines, 0.8) AS ptf,
         ST_EndPoint(lines) AS ptg
  FROM furthest_lines
), fl_dists AS ( -- Find distance from point to boundary
  SELECT fl_pts.id,
         ST_Distance(pta, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS da,
         ST_Distance(ptb, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS db,
         ST_Distance(ptc, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS dc,
         ST_Distance(ptd, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS dd,
         ST_Distance(pte, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS de,
         ST_Distance(ptf, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS df,
         ST_Distance(ptg, ST_ExteriorRing((ST_Dump(bd.geom)).geom)) AS dg
  FROM fl_pts, bd
  WHERE fl_pts.id = bd.id
), greatest AS ( -- Find greatest distance
  SELECT id, GREATEST(da, db, dc, dd, de, df, dg) AS radius
  FROM fl_dists
) -- Use this distance as radius for inscribed circle
UPDATE bd SET dm_inscr_c = 2 * radius
             FROM greatest
             WHERE bd.id = greatest.id;

DELETE FROM bd
  WHERE dm_inscr_c IS NULL;
```

### Compute the basic indicators
```sql
-- Perimeter
UPDATE bd SET b_perim = ROUND(CAST(ST_Perimeter(geom) AS NUMERIC), 2);

-- Area
UPDATE bd SET b_area = ROUND(CAST(ST_Area(geom) AS NUMERIC), 2);

-- Ajout volume
UPDATE bd SET b_r_vol_fac = b_area / b_perim;

-- Min and Max edges
WITH bd_segment AS (
  SELECT
    ST_PointN(geom, generate_series(1, ST_NPoints(geom)-1)) AS sp,
    ST_PointN(geom, generate_series(2, ST_NPoints(geom)  )) AS ep
  FROM
     -- extract the individual linestrings
     (SELECT (ST_Dump(ST_Boundary(geom))).geom
     FROM bd) AS linestrings
), bd_segment_geom AS (
  SELECT sp, ep, st_makeline(sp, ep) AS edge
  FROM bd_segment
), bd_segment_id AS (
  SELECT bd.id, ST_Length(ge.edge) AS length, ge.edge
  FROM bd_segment_geom ge
  JOIN bd ON ST_Touches(ge.edge, bd.geom)
  GROUP BY bd.id, ge.sp, ge.ep, ge.edge
), e_lgth AS (
  SELECT id, MAX(length) AS max, MIN(length) AS min FROM bd_segment_id
  GROUP BY id
)
UPDATE bd SET b_maxEdge = max,
              b_minEdge = min
          FROM e_lgth
          WHERE bd.id = e_lgth.id;

-- Number of stories
UPDATE bd SET b_stories = CASE
                            WHEN ROUND(b_height/3) < 1 THEN 1
                            ELSE ROUND(b_height/3)
                          END;

-- Floorspace
UPDATE bd SET b_floorsqm = b_area * b_stories;
```

### Compute the compacity indicators

```sql
-- Miller
UPDATE bd SET c_Miller = b_area / pow(.282 * b_perim, 2);
-- Schumm
UPDATE bd SET c_Schumm = 2 * sqrt(b_area / pi()) /  (2 * rad)
          FROM ccirc
          WHERE bd.id = ccirc.id;
-- Haggett TODO problem
UPDATE bd SET c_Haggett = dm_inscr_c / (2 * rad)
          FROM ccirc
          WHERE bd.id = ccirc.id;
-- Lee & Sallee
WITH cercleaire AS (
  SELECT bd.id, ST_Buffer(ST_Centroid(geom), sqrt(b_area / pi())) AS centcir FROM bd
), ops AS (
  SELECT bd.id, ST_Area(ST_Intersection(centcir, geom)) AS intersec,
         ST_Area(ST_Union(centcir, geom)) AS union_area
  FROM bd, cercleaire
  WHERE bd.id = cercleaire.id
)
UPDATE bd SET c_LeeSallee = intersec / union_area
       FROM ops
       WHERE bd.id = ops.id;
-- Ehrenburg
UPDATE bd SET c_Ehrenb = (pi() * pow(dm_inscr_c / 2, 2)) / b_area;
```

### Compute the bounding box indicators

```sql
-- Perim & area
UPDATE bd SET bb_perim = ST_Perimeter(bb),
              bb_area = ST_Area(bb)
          FROM bbox
          WHERE bd.id = bbox.id;
-- Width & length
WITH bboxcoord AS ( -- Extremes coords
 SELECT ST_XMin(bb) AS xmin,
        ST_XMax(bb) AS xmax,
        ST_YMin(bb) AS ymin,
        ST_YMax(bb) AS ymax,
        id
 FROM bbox
), matridist AS ( -- only 3 points needed
 SELECT ST_MakePoint(xmin, ymin) AS xymin,
        ST_MakePoint(xmin, ymax) AS xmiyma,
        ST_MakePoint(xmax, ymin) AS xmaymi,
        id
 FROM bboxcoord
), dist AS ( -- calculate both distances
 SELECT ST_Distance(xymin, xmiyma) AS distab,
        ST_Distance(xymin, xmaymi) AS distad,
        id
 FROM matridist
) -- Assign the longest as the length, shortest for the width
UPDATE bd SET bb_width = CASE
                           WHEN distab >= distad THEN distad
                           ELSE distab
                         END,
              bb_length = CASE
                            WHEN distab >= distad THEN distab
                            ELSE distad
                          END
          FROM dist
          WHERE bd.id = dist.id;
-- Ratios
UPDATE bd SET bb_r_lw = bb_length / bb_width,
              bb_r_area = b_area / bb_area,
              bb_r_perim = b_perim / bb_perim;
```

### Compute the circumscribed circle indicators

```sql
-- Radius
UPDATE bd SET cc_rad = rad
          FROM ccirc
          WHERE bd.id = ccirc.id;
-- Exchange index
WITH cercleaire AS ( -- Circle of same area centered on centroid
  SELECT id, ST_Buffer(ST_Centroid(geom), sqrt(b_area / pi())) AS centcir, geom
  FROM bd
), ops AS ( -- Area of intersection
  SELECT id, ST_Area(ST_Intersection(centcir, geom)) AS intersec
  FROM cercleaire
)
UPDATE bd SET cc_exch = intersec / b_area
             FROM ops
             WHERE bd.id = ops.id;
-- Detour index
UPDATE bd SET cc_detour = (2 * sqrt(pi() * b_area)) / ST_Perimeter(chull)
             FROM convhull
             WHERE bd.id = convhull.id;
```
### Compute the convex hull indicators
```sql
UPDATE bd SET ch_area = ST_Area(chull),
              ch_perim = ST_Perimeter(chull),
              ch_r_area = ST_Area(chull) / b_area,
              ch_r_perim = ST_Perimeter(chull) / b_perim
          FROM convhull
          WHERE bd.id = convhull.id;
```
### Compute the skeleton and centerline indicators

```sql
-- Number of deadends
-- Create table with skeleton's points
CREATE TABLE skelpt AS
 SELECT id, (ST_DumpPoints(skel)).geom
 FROM skeleton;
-- Add an id to the points
ALTER TABLE skelpt ADD COLUMN sid SERIAL PRIMARY KEY;

WITH s_unq_pts AS (-- Select unique points of skeleton
 SELECT DISTINCT id, ST_AsText(geom) FROM skelpt
), unq_pts_tot AS (-- Count them
 SELECT id, COUNT(*) AS unqpts
 FROM s_unq_pts
 GROUP BY id
), s_mlt_pts AS (-- Select points seen several times
 SELECT a.id, a.geom, a.sid
 FROM skelpt AS a, skelpt AS b
 WHERE ST_Equals(a.geom, b.geom) AND a.sid <> b.sid AND a.id = b.id
), s_mlt_diff AS (
 SELECT DISTINCT ST_AsText(geom), id
 FROM s_mlt_pts
), mlt_pts_tot AS (
 SELECT id, COUNT(*) AS mltpts
 FROM s_mlt_diff
 GROUP BY id
), deadends AS (
 SELECT unq_pts_tot.id AS id, unqpts - mltpts AS deadend
 FROM unq_pts_tot, mlt_pts_tot
 WHERE unq_pts_tot.id = mlt_pts_tot.id
)
UPDATE bd SET s_deadend = deadend
            FROM deadends
            WHERE bd.id = deadends.id;

-- Number of lines in centerline
WITH ctline AS (
  SELECT id, (ST_Dump(ctl)).geom
  FROM skeleton
), cttotline AS (
  SELECT id, COUNT(*) AS tot
  FROM ctline
  GROUP BY id
)
UPDATE bd SET sc_lines = tot
             FROM cttotline
             WHERE bd.id = cttotline.id;

-- Centerline length
UPDATE bd SET sc_length = ST_Length(ctl)
             FROM skeleton
             WHERE bd.id = skeleton.id;

 -- Number of centerline orientation
WITH nbor AS (
 SELECT id, COUNT(DISTINCT (id, orlc)) AS totor
 FROM sc_orien
 GROUP BY id
)
UPDATE bd SET sc_orient = totor
          FROM nbor
          WHERE bd.id = nbor.id;

-- Length of centerline by specific orientation
-- TODO must be a better way to update at once ?!
-- South - North
WITH sn AS (
  SELECT DISTINCT id, orlc, ltot
  FROM sum_orien
  WHERE orlc = 'SN'
)
UPDATE bd SET sc_l_sn = ltot
          FROM sn
          WHERE bd.id = sn.id;

-- East-West
WITH ew AS (
  SELECT DISTINCT id, orlc, ltot
  FROM sum_orien
  WHERE orlc = 'EW'
)
UPDATE bd SET sc_l_ew = ltot
          FROM ew
          WHERE bd.id = ew.id;

-- NorthEast - SouthWest
WITH nesw AS (
  SELECT DISTINCT id, orlc, ltot
  FROM sum_orien
  WHERE orlc = 'SWNE'
)
UPDATE bd SET sc_l_nesw = ltot
          FROM nesw
          WHERE bd.id = nesw.id;

-- SouthEast - NorthWest
WITH senw AS (
  SELECT DISTINCT id, orlc, ltot
  FROM sum_orien
  WHERE orlc = 'SENW'
)
UPDATE bd SET sc_l_senw = ltot
          FROM senw
          WHERE bd.id = senw.id;

-- Main orientation of centerline
WITH max_length AS (
  SELECT b.id, b.orlc, ltot
  FROM (
    SELECT id, orlc, MAX(ltot) OVER (PARTITION BY id) max_long
    FROM sum_orien
  ) a, sum_orien b
  WHERE ltot = max_long
), main_or AS (
  SELECT DISTINCT * FROM max_length
)
UPDATE bd SET sc_m_orient = orlc
          FROM main_or
          WHERE bd.id = main_or.id;
```

### Compute the miscellaneous indicators
```sql
-- Average distance to corners
WITH ptspoly AS (
  SELECT id, (ST_DumpPoints(geom)).geom AS ptsbld
  FROM bd
), corners AS ( -- Drop polygon closing point
  SELECT DISTINCT ptsbld AS corner, id  
  FROM ptspoly
), centblg AS (
  SELECT id, ST_Centroid(geom) AS center
  FROM bd
), dists AS (      
  SELECT corners.id, AVG(ST_Distance(corner, center)) AS dist
  FROM corners, centblg
  WHERE corners.id = centblg.id
  GROUP BY corners.id
)
UPDATE bd SET m_corndis = dist
          FROM dists
          WHERE bd.id = dists.id;

-- Number of inner courtyards
WITH holes AS (
  SELECT id, SUM(ST_NumInteriorRings(geom)) AS tot_holes
  FROM (SELECT id, (ST_Dump(geom)).geom As geom
  	FROM bd) AS a
  GROUP BY id
)
UPDATE bd SET m_court = tot_holes
          FROM holes
          WHERE bd.id = holes.id;

-- Relative area of courtyard
WITH rings AS (
  SELECT id,
         (ST_DumpRings((ST_Dump(geom)).geom)).path[1] as n,
         (ST_DumpRings((ST_Dump(geom)).geom)).geom as geom
  FROM bd
), holes AS (
  SELECT id, ST_Area(geom) AS area
  FROM rings
  WHERE n > 0
), holes_ar AS (
  SELECT id, SUM(area) AS tot_area
  FROM holes
  GROUP BY id
)
UPDATE bd SET m_court_area = tot_area
          FROM holes_ar
          WHERE bd.id = holes_ar.id;

UPDATE bd SET m_court_rel_a = m_court_area / b_area;


SELECT id, b_perim,b_area,b_r_vol_fac,b_MaxEdge,b_MinEdge,b_stories,b_floorsqm,c_Miller,c_Schumm,c_Haggett,c_LeeSallee,c_Ehrenb,bb_perim,bb_area,bb_length,bb_width,bb_r_lw,bb_r_area,bb_r_perim,cc_rad,cc_exch,cc_detour,ch_area,ch_perim,ch_r_area,ch_r_perim,s_deadend,sc_lines,sc_length,sc_orient,sc_l_sn,sc_l_ew,sc_l_nesw,sc_l_senw,m_corndis,m_court,m_court_area,m_court_rel_a, geom
INTO indic_bldg
FROM bd;
```
