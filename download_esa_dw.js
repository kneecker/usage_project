var translit = require('users/nikkuzmenko31/work:transliteration');
var russia_regions = ee.FeatureCollection(
  'projects/ee-nikkuzmenko31/assets/work_rus_regions_2023_with_buffer_1km');
var dw_palette = [
    '419BDF', '397D49', '88B053', '7A87C6', 'E49635',
    'DFC35A', 'C4281B', 'A59B8F', 'B39FE1'
];
var esa_palette =[
  '006400','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'ffbb22','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'ffff4c','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'f096ff','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'fa0000','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'b4b4b4','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  'f0f0f0','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  '0064c8','0000','0000','0000','0000','0000','0000','0000','0000','0000',
  '0096a0','0000','0000','0000','0000',
  '00cf75','0000','0000','0000','0000',
  'fae6a0'
];

var region_name = 'Краснодарский край';
var name_property = 'name';
var aoi = russia_regions.filter(ee.Filter.eq(name_property, region_name)).first().geometry();

var visualize = true;
var year = '2020';
var date_from = year + '-04-01';
var date_to = year + '-11-01';
var translit_name = translit.translit(region_name);
for(var i = 0; i < 10; i++){
    translit_name = translit_name.replace(' ', '_');
}

if (year == '2021' || year == '2020'){
  if (year == '2021'){
    var esa = ee.ImageCollection('ESA/WorldCover/v200').first();
    var esa_name = 'esa_2021';
  } else{
    var esa = ee.ImageCollection('ESA/WorldCover/v100').first();
    var esa_name = 'esa_2020';
  }
  var esa_crop = esa.clip(aoi);
    Export.image.toDrive({
    image: esa_crop,
    description: translit_name+'_'+esa_name,
    fileNamePrefix: translit_name+'_'+esa_name,
    scale: 10,
    crs: 'EPSG:3857',
    region: aoi,
    maxPixels: 1e13
  });
  if (visualize){
    Map.addLayer(
      esa_crop,
      {bands: ['Map'], min: 10, max: 100, palette: esa_palette},
      esa_name
    );
    Map.centerObject(aoi, 7);
  }
}
var year_int = parseInt(year, 10);
if (year_int > 2014){
  var dw_name = 'dw_' + year;
  var filter = ee.Filter.and(
      ee.Filter.bounds(aoi),
      ee.Filter.date(date_from, date_to));
  var dw = ee.ImageCollection('GOOGLE/DYNAMICWORLD/V1').filter(filter);
  var dw_crop = dw.select('label').mode().clip(aoi);
  Export.image.toDrive({
    image: dw_crop,
    description: translit_name+'_'+dw_name,
    fileNamePrefix: translit_name+'_'+dw_name,
    scale: 10,
    crs: 'EPSG:3857',
    region: aoi,
    maxPixels: 1e13
  });
  if (visualize){
    Map.addLayer(
      dw_crop,
      {min: 0, max: 8, palette: dw_palette, bands: ['label']},
      dw_name
    );
    Map.centerObject(aoi, 7);
  }
}
