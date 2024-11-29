extern crate retiles;

use core::num;

use convert_base::Convert;

fn base10_digits(number: u32) -> Vec<u8> {
    let mut digits = Vec::new();
    let mut n = number;
    while n > 0 {
        digits.push((n % 10) as u8);
        n /= 10;
    }
    digits
}

fn north_west_corner_same_integer_no_leading_0(number: u32) -> Vec<retiles::TileAssembly> {
    let mut ta = retiles::TileAssembly::new(retiles::TileSet::get_collatz_tileset());
    let base10_digits = base10_digits(number);
    let mut base2 = Convert::new(10, 2);
    let mut base3 = Convert::new(10, 3);

    let mut base2_digits = base2.convert::<u8, u8>(&base10_digits);
    let mut base3_digits = base3.convert::<u8, u8>(&base10_digits);
    base2_digits.reverse();
    base3_digits.reverse();

    for (i, base2_digit) in base2_digits.iter().enumerate() {
        //println!("{}: {}", i, base2_digit);
        ta.add_edge(&((i as i32, 0), retiles::Direction::RIGHT), *base2_digit);
    }

    for (i, base3_digit) in base3_digits.iter().enumerate() {
        //println!("{}: {}", i, base3_digit);
        ta.add_edge(
            &((0, -1 * i as i32), retiles::Direction::DOWN),
            *base3_digit,
        );
    }

    ta.solve_non_det(2)
}

use retiles::{Direction, Glue, NormalisedEdgePosition, RestrictedDirection, TileAssembly};

use svg::node::element::path::Data;
use svg::node::element::Path;
use svg::Document;

fn valued_edge_to_svg(edge: (&NormalisedEdgePosition, &Glue)) -> Vec<Path> {
    let scale = 2;
    fn restricted_direction_to_vec(d: RestrictedDirection, scale: i32) -> (i32, i32) {
        match d {
            RestrictedDirection::LEFT => (-1 * scale, 0),
            RestrictedDirection::UP => (0, -1 * scale),
        }
    }

    fn glue_to_color(g: &Glue) -> String {
        match g {
            0 => "#E6C58F".into(),
            1 => "#0AC52E".into(),
            2 => "#E37998".into(),
            _ => "#3FB3F3".into(),
        }
    }

    fn replace_0_with(t: (i32, i32), r: i32) -> (i32, i32) {
        match t {
            (0, 0) => (r, r),
            (0, _) => (r, t.1),
            (_, 0) => (t.0, r),
            t => (t.0, t.1),
        }
    }

    let starting_point = (scale * edge.0 .0 .0, -1 * scale * edge.0 .0 .1);
    let vec = restricted_direction_to_vec(edge.0 .1, scale);
    let replace_vec_1 = replace_0_with(vec, scale);
    let replace_vec_2 = replace_0_with(vec, -1 * scale);
    let half_vec_1 = (
        (replace_vec_1.0 as f32) * 0.5,
        (replace_vec_1.1 as f32) * 0.5,
    );
    let half_vec_2 = (
        (replace_vec_2.0 as f32) * 0.5,
        (replace_vec_2.1 as f32) * 0.5,
    );

    let svg_edge = Data::new()
        .move_to(starting_point.clone())
        .line_by(vec)
        .close();

    let svg_triangle_1 = Data::new()
        .move_to(starting_point)
        .line_by(half_vec_1)
        .line_by(half_vec_2)
        .close();

    let svg_triangle_2 = Data::new()
        .move_to(starting_point)
        .line_by(half_vec_2)
        .line_by(half_vec_1)
        .close();

    let stroke_width_multiplier = 0.02;

    let path_edge = Path::new()
        .set("stroke", "gray")
        .set("stroke-width", (scale as f32) * stroke_width_multiplier)
        .set("d", svg_edge);

    let path_triangle_1 = Path::new()
        .set("stroke", "gray")
        .set("fill", glue_to_color(edge.1))
        .set("stroke-width", (scale as f32) * stroke_width_multiplier)
        .set("d", svg_triangle_1);

    let path_triangle_2 = Path::new()
        .set("stroke", "gray")
        .set("fill", glue_to_color(edge.1))
        .set("stroke-width", (scale as f32) * stroke_width_multiplier)
        .set("d", svg_triangle_2);

    [path_edge, path_triangle_1, path_triangle_2].to_vec()
}

fn ta_to_svg(ta: &TileAssembly) -> String {
    let mut document = Document::new().set("viewBox", (-100, -100, 300, 300));

    let mut k = 0;
    for edge in ta.edges.iter() {
        let vec_path = valued_edge_to_svg(edge);

        for path in vec_path {
            document = document.add(path);
        }

        k += 1;
    }

    document.to_string()
}
use retiles::TileSet;
fn main() {
    let mut hydra_ta = TileAssembly::new(TileSet::get_collatz_tileset());
    hydra_ta.add_edge(&((0, 0), Direction::LEFT), 1);
    hydra_ta.add_edge(&((-1, 0), Direction::LEFT), 1);

    let N = 50;

    for i in 1..=N {
        hydra_ta.add_edge(&((-1 - i, 0), Direction::LEFT), 0);
    }

    let mut curr_position = (0, 0);

    let mut k = 0;
    while k <= N {
        eprintln!("k: {}", k);
        let curr_glue = hydra_ta
            .edges
            .get(&(curr_position, RestrictedDirection::LEFT))
            .unwrap();

        curr_position = (curr_position.0 - 1, curr_position.1 - 1);
        hydra_ta.add_edge(&(curr_position, Direction::UP), *curr_glue);
        hydra_ta.solve();

        k += 1;
    }

    println!("{}", ta_to_svg(&hydra_ta));
}
