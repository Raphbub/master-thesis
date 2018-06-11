# Data sets

The present work uses proprietary and open access data to perform the analyses.

Indicators calculated in the [analyse folder](https://github.com/Raphbub/master-thesis/tree/master/analyse) can be replicated using the sets presented below

## Swisstopo

#### TLM3D
The [TLM3D data set](https://shop.swisstopo.admin.ch/fr/products/landscape/tlm3D) has multiple types of features arranged in several categories. Three of those categories have been used in this analysis: the roads and ways, the buildings, and the lakes' contours. Samples of the data can be found [here](https://cms.geo.admin.ch/Topo/swisstlm3d/LV03/swissTLM3D_1.6_LV03_LN02_shp3d.zip)

#### swissBUILDINGS3D 1.0
The [swissBUILDINGS3D 1.0 data set](https://shop.swisstopo.admin.ch/fr/products/landscape/build3D) was used to attribute a height to the TLM3D buildings. It is no longer updated but samples can still be found [here](https://cms.geo.admin.ch/Topo/swissbuildings3d1/swissbuildings3dlv03.zip)

#### swissBOUNDARIES3D

The [swissBOUNDARIES3D data set](https://opendata.swiss/fr/dataset/swissboundaries3d-landesgrenzen) is freely available on the website for the swiss open data. It was used in the analysis to determine which municipalities were in the statistical agglomerations.

## OpenStreetMap
Data added in [OpenStreetMap](https://www.openstreetmap.org/#map=18/46.52590/6.58019) can be downloaded directly from the website or via third parties such as [Geofabrik](https://download.geofabrik.de/) or [OverPass Turbo](http://overpass-turbo.eu/). As the map is updated nightly, the more recent the download the more precise the data (technically...).
The OSM extract used for the analysis was downloaded on April 9 2018.

## Federal statistical office
Two data sets of the [national statistical office](https://www.bfs.admin.ch/bfs/en/home.html) were used to delineate the CCA agglomerations
- [Census statistics](https://www.bfs.admin.ch/bfs/fr/home/services/geostat/geodonnees-statistique-federale/btiments-logements-menages-personnes/resultats-recensement-depuis-2010.assetdetail.3543467.html) concerning the population at the hectare scale
- [Business census statistics](https://www.bfs.admin.ch/bfs/fr/home/services/geostat/geodonnees-statistique-federale/etablissements-emplois/statistique-structurel-entreprises-statent-depuis-2011.assetdetail.3303058.html) concerning the enterprises and available at the hectare level.

The final [data set](https://www.bfs.admin.ch/bfs/en/home/statistics/catalogues-databases/press-releases.assetdetail.188853.html) is about the urban character of space as defined in 2012 by the FSO. It was used to delineate the statistical agglomerations.
