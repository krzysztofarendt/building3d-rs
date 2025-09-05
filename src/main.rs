use building3d::Point;

fn main() {
    let pt = Point::new(1., 2.3, 4.56789);
    println!("{}", pt);
    println!("{:.3}", pt);
    println!("{:?}", pt);
}
