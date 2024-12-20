use std::{
    fmt::Display,
    ops::{Add, AddAssign, Index, IndexMut, Sub},
};

#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub struct Coord {
    pub x: i32,
    pub y: i32,
}

impl Coord {
    pub fn new(x: i32, y: i32) -> Self {
        Coord { x, y }
    }

    #[inline]
    pub fn rectilinear_distance(&self, other: &Self) -> i32 {
        (self - other).rectilinear_norm()
    }

    #[inline]
    pub fn rectilinear_norm(&self) -> i32 {
        self.x.abs() + self.y.abs()
    }

    #[inline]
    pub fn step_north(&self, distance: i32) -> Coord {
        Coord::new(self.x, self.y - distance)
    }
    #[inline]
    pub fn step_east(&self, distance: i32) -> Coord {
        Coord::new(self.x + distance, self.y)
    }
    #[inline]
    pub fn step_south(&self, distance: i32) -> Coord {
        Coord::new(self.x, self.y + distance)
    }
    #[inline]
    pub fn step_west(&self, distance: i32) -> Coord {
        Coord::new(self.x - distance, self.y)
    }

    #[inline]
    pub fn north(&self) -> Coord {
        Coord::new(self.x, self.y - 1)
    }
    #[inline]
    pub fn east(&self) -> Coord {
        Coord::new(self.x + 1, self.y)
    }

    #[inline]
    pub fn south(&self) -> Coord {
        Coord::new(self.x, self.y + 1)
    }

    #[inline]
    pub fn west(&self) -> Coord {
        Coord::new(self.x - 1, self.y)
    }

    #[inline]
    pub fn north_east(&self) -> Coord {
        Coord::new(self.x + 1, self.y - 1)
    }

    #[inline]
    pub fn north_west(&self) -> Coord {
        Coord::new(self.x - 1, self.y - 1)
    }

    #[inline]
    pub fn south_east(&self) -> Coord {
        Coord::new(self.x + 1, self.y + 1)
    }

    #[inline]
    pub fn south_west(&self) -> Coord {
        Coord::new(self.x - 1, self.y + 1)
    }

    #[inline]
    pub fn orthogs(&self) -> [Coord; 4] {
        [self.north(), self.east(), self.south(), self.west()]
    }
    #[inline]
    pub fn diags(&self) -> [Coord; 4] {
        [
            self.north_east(),
            self.south_east(),
            self.south_west(),
            self.north_west(),
        ]
    }
    #[inline]
    pub fn orthog_steps(&self, distance: i32) -> [Coord; 4] {
        [
            self.step_north(distance),
            self.step_east(distance),
            self.step_south(distance),
            self.step_west(distance),
        ]
    }
    #[inline]
    pub fn step_return(&mut self, other: &Self) -> Self {
        let prev = *self;
        *self += *other;
        prev
    }
}

impl Add for Coord {
    type Output = Coord;
    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        Coord {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Add for &Coord {
    type Output = Coord;

    #[inline]
    fn add(self, rhs: Self) -> Self::Output {
        *self + *rhs
    }
}

impl AddAssign for Coord {
    #[inline]
    fn add_assign(&mut self, rhs: Self) {
        self.x += rhs.x;
        self.y += rhs.y;
    }
}

impl Sub for &Coord {
    type Output = Coord;
    #[inline]
    fn sub(self, rhs: Self) -> Self::Output {
        Coord::new(self.x - rhs.x, self.y - rhs.y)
    }
}

#[derive(Debug, Clone)]
pub struct Grid<T> {
    pub items: Vec<T>,
    pub width: usize,
    pub height: usize,
}

impl<T> Index<Coord> for Grid<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: Coord) -> &Self::Output {
        &self.items[self.width * index.y as usize + index.x as usize]
    }
}

impl<T> IndexMut<Coord> for Grid<T> {
    #[inline]
    fn index_mut(&mut self, index: Coord) -> &mut Self::Output {
        &mut self.items[self.width * index.y as usize + index.x as usize]
    }
}

impl<T: Copy> Grid<T> {
    pub fn new(item: T, width: usize, height: usize) -> Grid<T> {
        Grid::from_iter(&mut (0..(width * height)).map(|_| item), width)
    }
}

impl<T: Display> Grid<T> {
    pub fn print(&self) {
        for (c, i) in self.iter() {
            print!("{i}");
            if c.x == self.width as i32 - 1 {
                println!();
            }
        }
    }

    pub fn print_func(&self, func: fn(&T) -> String) {
        for (c, i) in self.iter() {
            print!("{}", func(i));
            if c.x == self.width as i32 - 1 {
                println!();
            }
        }
    }

    pub fn print_with_commas(&self) {
        for (c, i) in self.iter() {
            print!("{i}");
            if c.x == self.width as i32 - 1 {
                println!();
            } else {
                print!(",");
            }
        }
    }
}
impl<T> Grid<T> {
    pub fn from_iter(iter: &mut dyn Iterator<Item = T>, width: usize) -> Self {
        let items: Vec<T> = iter.collect();
        let height = items.len() / width;
        Grid {
            items,
            width,
            height,
        }
    }

    #[inline]
    pub fn bounds_check(&self, coord: &Coord) -> bool {
        coord.x < self.width as i32 && coord.y < self.height as i32 && coord.x >= 0 && coord.y >= 0
    }

    pub fn iter<'a>(
        &'a self,
    ) -> std::iter::Map<
        std::iter::Enumerate<std::slice::Iter<'a, T>>,
        impl FnMut((usize, &'a T)) -> (Coord, &T),
    > {
        let width = self.width;
        self.items
            .iter()
            .enumerate()
            .map(move |(i, item)| (Coord::new((i % width) as i32, (i / width) as i32), item))
    }

    pub fn coord_iter(&self) -> std::iter::Map<std::ops::Range<usize>, impl FnMut(usize) -> Coord> {
        let width = self.width;
        (0..self.items.len()).map(move |i| Coord::new((i % width) as i32, (i / width) as i32))
    }

    #[inline]
    pub fn orthogs_coords(&self, coord: &Coord) -> [Option<(Coord, &T)>; 4] {
        coord
            .orthogs()
            .map(|c| self.checked_index(&c).map(|item| (c, item)))
    }

    #[inline]
    pub fn diags(&self, coord: &Coord) -> [Option<&T>; 4] {
        coord.diags().map(|c| self.checked_index(&c))
    }

    #[inline]
    pub fn diags_unchecked(&self, coord: &Coord) -> [&T; 4] {
        coord.diags().map(|c| self.index(c))
    }

    #[inline]
    pub fn checked_index(&self, coord: &Coord) -> Option<&T> {
        if self.bounds_check(coord) {
            Some(&self[*coord])
        } else {
            None
        }
    }
}

impl<T: PartialEq + Copy> Grid<T> {
    pub fn iter_filtered(
        &self,
        value: T,
    ) -> std::iter::Map<
        std::iter::Filter<
            std::iter::Enumerate<std::iter::Copied<std::slice::Iter<'_, T>>>,
            impl FnMut(&(usize, T)) -> bool,
        >,
        impl FnMut((usize, T)) -> (Coord, T),
    > {
        let width = self.width;
        self.items
            .iter()
            .copied()
            .enumerate()
            .filter(move |(_, item)| *item == value)
            .map(move |(i, item)| (Coord::new((i % width) as i32, (i / width) as i32), item))
    }
}

impl<'a, T> Grid<T> {
    pub fn stride_iter(&'a self, start: Coord, stride: Coord) -> GridLineIter<'a, T> {
        GridLineIter {
            grid: self,
            coord: start,
            stride,
        }
    }

    pub fn stride_iter_mut(&'a mut self, start: Coord, stride: Coord) -> GridLineIterMut<'a, T> {
        GridLineIterMut {
            grid: self,
            coord: start,
            stride,
        }
    }

    pub fn north_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(0, -1))
    }

    pub fn east_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(1, 0))
    }

    pub fn south_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(0, 1))
    }

    pub fn west_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(-1, 0))
    }

    pub fn north_east_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(1, -1))
    }

    pub fn north_west_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(-1, -1))
    }

    pub fn south_east_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(1, 1))
    }

    pub fn south_west_iter(&'a self, start: Coord) -> GridLineIter<'a, T> {
        self.stride_iter(start, Coord::new(-1, 1))
    }
    pub fn orthogonal_direction_iterators(&'a self, start: Coord) -> [GridLineIter<'a, T>; 4] {
        [
            self.north_iter(start),
            self.east_iter(start),
            self.south_iter(start),
            self.west_iter(start),
        ]
    }
    pub fn diagonal_direction_iterators(&'a self, start: Coord) -> [GridLineIter<'a, T>; 4] {
        [
            self.north_east_iter(start),
            self.south_east_iter(start),
            self.south_west_iter(start),
            self.north_west_iter(start),
        ]
    }
    pub fn all_direction_iterators(&'a self, start: Coord) -> [GridLineIter<'a, T>; 8] {
        [
            self.north_iter(start),
            self.east_iter(start),
            self.south_iter(start),
            self.west_iter(start),
            self.north_east_iter(start),
            self.south_east_iter(start),
            self.south_west_iter(start),
            self.north_west_iter(start),
        ]
    }

    pub fn north_iter_mut(&'a mut self, start: Coord) -> GridLineIterMut<'a, T> {
        self.stride_iter_mut(start, Coord::new(0, -1))
    }

    pub fn east_iter_mut(&'a mut self, start: Coord) -> GridLineIterMut<'a, T> {
        self.stride_iter_mut(start, Coord::new(1, 0))
    }

    pub fn south_iter_mut(&'a mut self, start: Coord) -> GridLineIterMut<'a, T> {
        self.stride_iter_mut(start, Coord::new(0, 1))
    }

    pub fn west_iter_mut(&'a mut self, start: Coord) -> GridLineIterMut<'a, T> {
        self.stride_iter_mut(start, Coord::new(-1, 0))
    }
}

pub struct GridLineIter<'a, T> {
    grid: &'a Grid<T>,
    coord: Coord,
    stride: Coord,
}

impl<'a, T> Iterator for GridLineIter<'a, T> {
    type Item = (&'a T, Coord);

    fn next(&mut self) -> Option<Self::Item> {
        if self.grid.bounds_check(&self.coord) {
            let coord = self.coord.step_return(&self.stride);
            Some((&self.grid[coord], coord))
        } else {
            None
        }
    }
}

pub struct GridLineIterMut<'a, T> {
    grid: &'a mut Grid<T>,
    coord: Coord,
    stride: Coord,
}

impl<'a, T> Iterator for GridLineIterMut<'a, T> {
    type Item = (&'a mut T, Coord);

    fn next(&mut self) -> Option<Self::Item> {
        if self.grid.bounds_check(&self.coord) {
            let coord = self.coord.step_return(&self.stride);
            let ptr = (&mut self.grid[coord]) as *mut T;
            unsafe { Some((&mut *ptr, coord)) }
        } else {
            None
        }
    }
}

#[cfg(test)]
mod test {
    use super::Grid;

    #[test]
    fn grid_iter() {
        let grid = Grid::from_iter(&mut (0..100u32), 10);
        for (c, item) in grid.iter() {
            println!("{c:?}, {item}");
        }
    }

    #[test]
    fn grid_iter_filtered() {
        let grid = Grid::from_iter(&mut (0..100u32), 10);
        for (c, item) in grid.iter_filtered(5) {
            println!("{c:?}, {item}");
        }
    }
}
