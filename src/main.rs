use serde::Serialize;
use serde::{Deserialize };
use std::fs;
use geo::{BooleanOps, BoundingRect, GeodesicArea, Intersects};
use geo::geometry::{Polygon, LineString};
use rstar::{RTree, AABB};
use std::time::Instant;


fn create_spatial_index(polygons: &Vec<Polygon>) -> RTree<Polygon> {
	RTree::bulk_load(polygons.clone())
}

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
	pub overlap_type: OverlapType
}

#[derive(Clone, Debug, Deserialize, Serialize)]
pub enum OverlapType {
	NoOverlap,
	MinorOverlap,
	MajorOverlap,
	DuplicateGeometry,
}

fn main() {
	let start_time = Instant::now();
	let polygons = get_coords_from_geojson("overlaps.geojson");
	println!("Number of polygons: {}", polygons.len());

	let irs = is_intersected_polygon(&polygons);
	for ir in irs {
		println!(
			"Intersection: Index: {}, \nintersected_with: {}, \nIntersects: {}, \nPercentage: {:?}%, \nOriginal Area: {:.4} ha, \nIntersected Area: {:.4} ha. \noverlap_type: {:?} \n\n",
			ir.index,
			ir.intersect_with,
			ir.is_intersect,
			ir.percentage_area,
			ir.original_area,
			ir.intersected_area,
			ir.overlap_type
		);
	}
	let end_time = start_time.elapsed();
	println!("Time elapsed is: {:?}", end_time);
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


pub fn classify_overlap(percentage: f64) -> OverlapType {
	match percentage {
		0.0 => OverlapType::NoOverlap,
		p if p > 0.00 && p < 10.00 => OverlapType::MinorOverlap,
		p if p >= 10.00 && p < 100.00 => OverlapType::MajorOverlap,
		p if p >= 99.99 => OverlapType::DuplicateGeometry,
		_ => OverlapType::NoOverlap,
	}
}
pub fn is_intersected_polygon(polygons: &Vec<Polygon>) -> Vec<IntersectionResult> {
	let mut intersection_result = Vec::new();

	let s_index = create_spatial_index(polygons);

	for (idx, i) in polygons.iter().enumerate() {
		let bbox = i.bounding_rect().unwrap();

		let candidates = s_index.locate_in_envelope(&AABB::from_corners(
			[bbox.min().x, bbox.min().y].into(),
			[bbox.max().x, bbox.max().y].into()
		));

		for (idxj, j) in candidates.enumerate() {
			if idx == idxj {
				continue;
			}

			let intersected = i.intersects(j);

			if intersected {
				let intersection = i.intersection(j);
				let intersection_area = intersection.geodesic_area_signed().abs() / 10_000.00;
				let is_intersected = intersection_area > 0.00;

				let polygon_a_area = i.geodesic_area_signed().abs() / 10_000.00;
				let intersected_a_percentage =( intersection_area/polygon_a_area) * 100.00;

				let results = IntersectionResult{
					index: idx,
					intersect_with: idxj,
					is_intersect: is_intersected,
					percentage_area: intersected_a_percentage,
					original_area: polygon_a_area,
					intersected_area: intersection_area,
					overlap_type: classify_overlap(intersected_a_percentage)
				};

				intersection_result.push(results)
			}
		}
	}
	intersection_result
}

pub fn float_precision(float: f64, precision: i32) -> f64 {
	let factor = 10_f64.powi(precision);
	(factor * float).round() / float
}