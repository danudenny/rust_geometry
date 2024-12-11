use serde::Serialize;
use serde::{Deserialize };
use std::fs;
use geo::{BooleanOps, GeodesicArea, Intersects};
use geo::geometry::{Polygon, LineString};

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Geometry {
	pub r#type: Option<String>,
	pub coordinates: Option<Vec<Vec<Vec<f64>>>>,
}
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Features {
	pub r#type: Option<String>,
	pub geometry: Option<Geometry>,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Root {
	pub r#type: Option<String>,
	pub features: Option<Vec<Features>>,
}


#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct IntersectionResult {
	pub index: usize,
	pub intersect_with: usize,
	pub is_intersect: bool,
	pub percentage_area: f64,
	pub original_area: f64,
	pub intersected_area: f64,
}

fn main() {
	let polygons = get_coords_from_geojson("overlaps.geojson");
	println!("Number of polygons: {}", polygons.len());
	// let area = get_area_sqm(&polygons);
	// println!("Area sqm: {:?}", area);
	// let area_ha = get_area_ha(&polygons);
	// println!("Area ha: {:?}", area_ha);
	let irs = is_intersected_polygon(&polygons);
	for ir in irs {
		println!(
			"Intersection: Index: {}, \nintersected_with: {}, \nIntersects: {}, \nPercentage: {:.2}%, \nOriginal Area: {:.4} ha, \nIntersected Area: {:.4} ha. \n\n",
			ir.index,
			ir.intersect_with,
			ir.is_intersect,
			ir.percentage_area,
			ir.original_area,
			ir.intersected_area,
		);
	}
}

pub fn get_coords_from_geojson(filename: &str) -> Vec<Polygon> {
	let file = fs::read_to_string(filename).unwrap();
	let geojson = serde_json::from_str::<Root>(&file).unwrap();

	let mut results = Vec::new();
	if let Some(features) = geojson.features {
		for feature in features {
			let geometry = feature.geometry.unwrap();
			let coords = geometry.coordinates.unwrap()[0].clone();
			let coords_tupl: Vec<(f64, f64)> = coords
				.into_iter()
				.map(|c| (c[0], c[1]))
				.collect();

			let polygon = Polygon::new(LineString::from(coords_tupl), vec![]);

			results.push(polygon.clone());
		}
	}
	results
}

pub fn get_area_sqm(polygons: &Vec<Polygon>) -> Vec<f64>{
	let mut areas = Vec::new();
	for polygon in polygons {
		areas.push(polygon.geodesic_area_signed().abs())
	}
	areas
}

pub fn get_area_ha(polygons: &Vec<Polygon>) -> Vec<f64>{
	let mut areas = Vec::new();
	for polygon in polygons {
		areas.push(polygon.geodesic_area_signed().abs() / 10_000_f64)
	}
	areas
}

pub fn is_intersected_polygon(polygons: &Vec<Polygon>) -> Vec<IntersectionResult> {
	let mut intersection_result = Vec::new();

	for (idx, i) in polygons.iter().enumerate() {
		for (idxj, j) in polygons.iter().enumerate() {
			if idx == idxj {
				continue;
			}

			let intersected = i.intersects(j);

			if intersected {
				let intersection = i.intersection(j);
				let intersection_area = intersection.geodesic_area_signed().abs() / 10_000_f64;
				let is_intersected = intersection_area > 0_f64;

				let polygon_a_area = i.geodesic_area_signed().abs() / 10_000_f64;
				let intersected_a_percentage =( intersection_area/polygon_a_area) * 100_f64;

				let results = IntersectionResult{
					index: idx,
					intersect_with: idxj,
					is_intersect: is_intersected,
					percentage_area: intersected_a_percentage,
					original_area: polygon_a_area,
					intersected_area: intersection_area
				};

				intersection_result.push(results)
			}
		}
	}
	intersection_result
}

