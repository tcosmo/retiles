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

fn main() {
    let number = 42;
    let tas = north_west_corner_same_integer_no_leading_0(number);
    println!("{}: {}", number, tas.len());
    println!("{}", serde_json::to_string(&tas[10]).unwrap());

    // for number in 1..100 {
    //     let tas = north_west_corner_same_integer_no_leading_0(number);
    //     println!("{}: {}", number, tas.len());
    // }
}
