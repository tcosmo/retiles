use serde::{Deserialize, Serialize};
use serde_with::serde_as;
use std::collections::{HashMap, HashSet, VecDeque};
use vector2math::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum Direction {
    UP,
    RIGHT,
    DOWN,
    LEFT,
}

/// A restricted direction is either UP or LEFT. Useful for NormalisedEdgePosition.
#[derive(Serialize, Deserialize, Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum RestrictedDirection {
    UP,
    LEFT,
}

const UP: (i32, i32) = (0, 1);
const RIGHT: (i32, i32) = (1, 0);
const DOWN: (i32, i32) = (0, -1);
const LEFT: (i32, i32) = (-1, 0);

pub type Glue = u8;

#[derive(Serialize, Deserialize, Copy, Clone, Debug)]
pub struct TileType {
    up: Glue,
    right: Glue,
    down: Glue,
    left: Glue,
}

impl TileType {
    /// Returns true if the tile type matches the given adjacent edges.
    pub fn is_matching(self, adjacent_edges: &AdjacentEdges) -> bool {
        if let Some(glue) = adjacent_edges.up.1 {
            if glue != self.up {
                return false;
            }
        }
        if let Some(glue) = adjacent_edges.right.1 {
            if glue != self.right {
                return false;
            }
        }
        if let Some(glue) = adjacent_edges.down.1 {
            if glue != self.down {
                return false;
            }
        }
        if let Some(glue) = adjacent_edges.left.1 {
            if glue != self.left {
                return false;
            }
        }
        true
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TileSet {
    tile_types: Vec<TileType>,
}

impl TileSet {
    /// Returns the tile types that match the given adjacent edges.
    pub fn matching_tiles(&self, adjacent_edges: &AdjacentEdges) -> Vec<usize> {
        let mut to_return: Vec<usize> = Vec::new();
        for (i, tile_type) in self.tile_types.iter().enumerate() {
            if tile_type.is_matching(adjacent_edges) {
                to_return.push(i);
            }
        }
        to_return
    }
    /// Returns the Collatz tileset, cf. Tristan's thesis
    /// # Examples
    ///  
    /// ```
    /// use retiles::get_collatz_tileset;
    /// let tileset = get_collatz_tileset();
    /// ```
    pub fn get_collatz_tileset() -> Self {
        TileSet {
            tile_types: vec![
                TileType {
                    up: 0,
                    right: 0,
                    down: 0,
                    left: 0,
                },
                TileType {
                    up: 0,
                    right: 1,
                    down: 1,
                    left: 0,
                },
                TileType {
                    up: 0,
                    right: 2,
                    down: 0,
                    left: 1,
                },
                TileType {
                    up: 1,
                    right: 0,
                    down: 1,
                    left: 1,
                },
                TileType {
                    up: 1,
                    right: 1,
                    down: 0,
                    left: 2,
                },
                TileType {
                    up: 1,
                    right: 2,
                    down: 1,
                    left: 2,
                },
            ],
        }
    }
}

pub type TilePosition = (i32, i32);
/// We denote edge by pairs such as ((0,0), LEFT)
/// The relation between TilePosition and EdgePosition is that
/// the right edge of the tile at TilePosition(0,0) is ((0,0), (0,-1)) = EdgePosition((0,0), LEFT).
pub type EdgePosition = ((i32, i32), Direction);
// Normalised edge position either has direction LEFT or UP.
pub type NormalisedEdgePosition = ((i32, i32), RestrictedDirection);

pub fn unormalise_edge_position(normalised_edge_position: &NormalisedEdgePosition) -> EdgePosition {
    let (base_pos, restricted_direction) = *normalised_edge_position;
    match restricted_direction {
        RestrictedDirection::UP => (base_pos, Direction::UP),
        RestrictedDirection::LEFT => (base_pos, Direction::LEFT),
    }
}

/// The edge (0,0) -> (-1,0) can be defined either by ((0,0), LEFT) or ((-1,0), RIGHT).
/// The normalised edge position is the one where the direction is always LEFT or UP.
pub fn normalise_edge_position(edge_position: &EdgePosition) -> NormalisedEdgePosition {
    let (base_pos, direction) = *edge_position;
    match direction {
        Direction::UP => (base_pos, RestrictedDirection::UP),
        Direction::LEFT => (base_pos, RestrictedDirection::LEFT),
        Direction::RIGHT => (base_pos.add(RIGHT), RestrictedDirection::LEFT),
        Direction::DOWN => (base_pos.add(DOWN), RestrictedDirection::UP),
    }
}
pub struct AdjacentEdges {
    up: (NormalisedEdgePosition, Option<Glue>),
    right: (NormalisedEdgePosition, Option<Glue>),
    down: (NormalisedEdgePosition, Option<Glue>),
    left: (NormalisedEdgePosition, Option<Glue>),
}

impl AdjacentEdges {
    pub fn number_of_defined_edges(&self) -> u8 {
        let mut to_return = 0;
        if self.up.1.is_some() {
            to_return += 1;
        }
        if self.right.1.is_some() {
            to_return += 1;
        }
        if self.down.1.is_some() {
            to_return += 1;
        }
        if self.left.1.is_some() {
            to_return += 1;
        }
        to_return
    }
}

/// Returns the tile positions that are adjacent to the given edge position.
pub fn adjacent_tile_positions(edge_position: NormalisedEdgePosition) -> [TilePosition; 2] {
    match edge_position {
        (base_pos, RestrictedDirection::UP) => [base_pos, base_pos.add(RIGHT)],
        (base_pos, RestrictedDirection::LEFT) => [base_pos, base_pos.add(DOWN)],
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum TileTypeOrImpossible {
    TileType(TileType),
    NoTileTypeCanFit,
}

#[serde_as]
#[derive(Serialize, Debug, Clone)]
/// A tile assembly is a set of tiles positioned on a grid.
pub struct TileAssembly {
    #[serde_as(as = "Vec<(_, _)>")]
    pub tiles: HashMap<TilePosition, TileTypeOrImpossible>,
    #[serde_as(as = "Vec<(_, _)>")]
    pub edges: HashMap<NormalisedEdgePosition, Glue>,
    pub tileset: TileSet,
    /// The frontier is the set of tile positions that have at least one edge defined.
    current_frontier: HashSet<TilePosition>,
}

impl TileAssembly {
    /// Returns all positions where more than 1 tile can fit.
    ///
    /// # Arguments
    /// * `threshold` - The number of matching edges that a tile type must have to be considered to fit.
    ///
    pub fn get_nondet_position(&self, threshold: u8) -> Vec<TilePosition> {
        let mut to_return: Vec<TilePosition> = Vec::new();
        for tile_position in self.current_frontier.iter() {
            let adjacent_edges = self.adjacent_edges(tile_position);

            if adjacent_edges.number_of_defined_edges() < threshold {
                continue;
            }

            let matching_tiles = self.tileset.matching_tiles(&adjacent_edges);
            if matching_tiles.len() > 1 {
                to_return.push(*tile_position);
            }
        }
        to_return
    }

    /// Similar to solve but returns all terminal assemblies.
    pub fn solve_non_det(&self, threshold: u8) -> Vec<TileAssembly> {
        let mut to_visit: VecDeque<TileAssembly> = VecDeque::new();
        to_visit.push_back(self.clone());
        let mut terminal_assemblies: Vec<TileAssembly> = Vec::new();

        while !to_visit.is_empty() {
            let mut current_assembly = to_visit.pop_front().unwrap();

            let old_frontier = current_assembly.current_frontier.clone();

            //println!("Solving frontier");

            current_assembly.solve_frontier();

            let non_det_positions = current_assembly.get_nondet_position(threshold);

            if old_frontier == current_assembly.current_frontier && non_det_positions.is_empty() {
                terminal_assemblies.push(current_assembly);
                continue;
            }

            if non_det_positions.is_empty() {
                to_visit.push_back(current_assembly);
                continue;
            }

            let non_det_position = non_det_positions[0];
            let adjacent_edges = current_assembly.adjacent_edges(&non_det_position);
            let matching_tiles = current_assembly.tileset.matching_tiles(&adjacent_edges);
            for tile_type_index in matching_tiles {
                let mut new_assembly = current_assembly.clone();
                // println!(
                //     "Non det: Adding tile {:?} at position {:?}",
                //     tile_type_index, non_det_position
                // );
                new_assembly
                    .add_tile_from_tileset(&non_det_position, tile_type_index)
                    .unwrap();
                new_assembly.current_frontier.remove(&non_det_position);
                to_visit.push_back(new_assembly);
            }
        }

        terminal_assemblies
    }

    /// Solves the tile assembly deterministically, i.e. positions where more than 1 tile can fit are left empty.
    pub fn solve(&mut self) {
        loop {
            let old_frontier = self.current_frontier.clone();
            self.solve_frontier();
            if old_frontier == self.current_frontier {
                break;
            }
        }
    }

    /// Solves the frontier of the tile assembly.
    pub fn solve_frontier(&mut self) {
        let mut positions_to_remove: Vec<TilePosition> = Vec::new();
        for tile_position in self.current_frontier.clone() {
            let adjacent_edges = self.adjacent_edges(&tile_position);
            let matching_tiles = self.tileset.matching_tiles(&adjacent_edges);
            if matching_tiles.is_empty() {
                //println!("Adding impossible tile at position {:?}", tile_position);
                self.tiles
                    .insert(tile_position, TileTypeOrImpossible::NoTileTypeCanFit);
            } else if matching_tiles.len() == 1 {
                let tile_type_index = matching_tiles[0];
                // println!(
                //     "Adding tile {:?} at position {:?}",
                //     tile_type_index, tile_position
                // );
                self.add_tile_from_tileset(&tile_position, tile_type_index)
                    .unwrap();
                positions_to_remove.push(tile_position);
            }
        }
        for position in positions_to_remove {
            self.current_frontier.remove(&position);
        }
    }

    /// Returns the edges that are adjacent to the given tile position in order UP, RIGHT, DOWN, LEFT.
    pub fn adjacent_edges(&self, tile_position: &TilePosition) -> AdjacentEdges {
        let up_edge_position = (tile_position.add(UP), RestrictedDirection::LEFT);
        let right_edge_position = (*tile_position, RestrictedDirection::UP);
        let down_edge_position = (*tile_position, RestrictedDirection::LEFT);
        let left_edge_position = (tile_position.add(LEFT), RestrictedDirection::UP);
        AdjacentEdges {
            up: (up_edge_position, self.edges.get(&up_edge_position).copied()),
            right: (
                right_edge_position,
                self.edges.get(&right_edge_position).copied(),
            ),
            down: (
                down_edge_position,
                self.edges.get(&down_edge_position).copied(),
            ),
            left: (
                left_edge_position,
                self.edges.get(&left_edge_position).copied(),
            ),
        }
    }

    pub fn add_normalised_edge(
        &mut self,
        edge_position: &NormalisedEdgePosition,
        glue: Glue,
    ) -> Result<(), ()> {
        self.add_edge(&unormalise_edge_position(edge_position), glue)
    }
    /// Adds an edge to the tile assembly.
    ///
    /// # Returns
    /// - Err(()) if the edge cannot be added (e.g. because it is already there)
    /// - OK(()) if the edge was added successfully
    pub fn add_edge(&mut self, edge_position: &EdgePosition, glue: Glue) -> Result<(), ()> {
        let edge_position = normalise_edge_position(edge_position);
        if let Some(other_glue) = self.edges.get(&edge_position) {
            if *other_glue != glue {
                return Result::Err(());
            }
        }

        self.edges.insert(edge_position, glue);
        for tile_position in adjacent_tile_positions(edge_position) {
            match self.tiles.get(&tile_position) {
                Some(TileTypeOrImpossible::NoTileTypeCanFit) => {
                    // This case where a side of the edge cannot fit any tile should be rare in practice but we allow it.
                    self.current_frontier.insert(tile_position);
                }
                None => {
                    self.current_frontier.insert(tile_position);
                }
                _ => {}
            }
        }
        Result::Ok(())
    }

    /// Returns true if the tile can be placed at the given position.
    /// This means that there is no tile at the position and that all
    /// potentially already defined edges match.
    ///
    pub fn is_tile_placeable_at_position(
        &self,
        tile_position: &TilePosition,
        tile_type: &TileType,
    ) -> bool {
        if self.tiles.contains_key(tile_position) {
            //println!("Tile already there");
            return false;
        }
        let adjacent_edges = self.adjacent_edges(tile_position);
        if let Some(&glue) = self.edges.get(&adjacent_edges.up.0) {
            if glue != tile_type.up {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.right.0) {
            if glue != tile_type.right {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.down.0) {
            if glue != tile_type.down {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.left.0) {
            if glue != tile_type.left {
                return false;
            }
        }
        true
    }

    /// Add a tile to the tile assembly, selected from the tile set.
    ///
    /// # Returns
    /// - Err(()) if the tile cannot be added (e.g. because it is already there or one of the edges that is already defined mismatches)
    /// - OK(()) if the tile was added successfully
    pub fn add_tile_from_tileset(
        &mut self,
        tile_position: &TilePosition,
        tile_type_index: usize,
    ) -> Result<(), ()> {
        let tile_type = self.tileset.tile_types[tile_type_index];
        if !self.is_tile_placeable_at_position(tile_position, &tile_type) {
            return Result::Err(());
        }

        let adjacent_edges = self.adjacent_edges(tile_position);

        let _ = self.add_normalised_edge(&adjacent_edges.up.0, tile_type.up);
        let _ = self.add_normalised_edge(&adjacent_edges.right.0, tile_type.right);
        let _ = self.add_normalised_edge(&adjacent_edges.down.0, tile_type.down);
        let _ = self.add_normalised_edge(&adjacent_edges.left.0, tile_type.left);

        self.tiles
            .insert(*tile_position, TileTypeOrImpossible::TileType(tile_type));

        Result::Ok(())
    }

    pub fn new(tileset: TileSet) -> Self {
        TileAssembly {
            tiles: HashMap::new(),
            edges: HashMap::new(),
            tileset,
            current_frontier: HashSet::new(),
        }
    }

    pub fn assembly_from_Collatz_parity_vector(parity_vector: Vec<u8>) -> Self {
        let mut to_return: TileAssembly = TileAssembly::new(TileSet::get_collatz_tileset());
        let mut current_position: TilePosition = (0, 0);
        for parity in parity_vector {
            if parity == 0 {
                let edge_to_add = (current_position, Direction::LEFT);
                //println!("Adding edge {:?}", edge_to_add);
                to_return.add_edge(&edge_to_add, 0).unwrap();
                current_position = current_position.add(LEFT);
            } else {
                let edge_to_add = (current_position, Direction::DOWN);
                //println!("Adding edge {:?}", edge_to_add);
                to_return.add_edge(&edge_to_add, 1).unwrap();
                current_position = current_position.add(DOWN);

                let edge_to_add = (current_position, Direction::LEFT);
                //println!("Adding edge {:?}", edge_to_add);
                to_return.add_edge(&edge_to_add, 0).unwrap();
                current_position = current_position.add(LEFT);
            }
        }
        to_return
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solve_non_det() {
        let mut ta = TileAssembly::new(TileSet::get_collatz_tileset());
        let _ = ta.add_edge(&((0, 0), Direction::LEFT), 0);
        let _ = ta.add_edge(&((-1, 0), Direction::DOWN), 0);

        for (tile_position, tile) in ta.tiles.iter() {
            println!("{:?} {:?}", tile_position, tile);
        }

        println!("{:?}", ta.solve_non_det(2).len());

        let mut ta = TileAssembly::new(TileSet::get_collatz_tileset());
        let _ = ta.add_edge(&((0, 0), Direction::RIGHT), 1);
        let _ = ta.add_edge(&((0, 1), Direction::RIGHT), 0);
        let _ = ta.add_edge(&((0, 2), Direction::RIGHT), 0);
        let _ = ta.add_edge(&((0, 3), Direction::RIGHT), 1);
        let _ = ta.add_edge(&((0, 4), Direction::RIGHT), 1);

        let _ = ta.add_edge(&((0, 0), Direction::DOWN), 2);
        let _ = ta.add_edge(&((0, -1), Direction::DOWN), 0);
        let _ = ta.add_edge(&((0, -2), Direction::DOWN), 1);

        for (tile_position, tile) in ta.tiles.iter() {
            println!("{:?} {:?}", tile_position, tile);
        }

        println!("{:?}", ta.solve_non_det(2).len());

        assert!(true);
    }

    #[test]
    fn initial_frontier_from_collatz_parity_vector() {
        let mut ta = TileAssembly::assembly_from_Collatz_parity_vector(vec![0, 1, 0, 1, 1]);

        for (edge_position, glue) in ta.edges.iter() {
            println!("{:?} {:?}", edge_position, glue);
        }

        for position in ta.current_frontier.iter() {
            println!("{:?}", position);
        }

        ta.solve();

        for (tile_position, tile_type_or_impossible) in ta.tiles.iter() {
            println!("{:?} {:?}", tile_position, tile_type_or_impossible);
        }

        assert_eq!(1, 1);
    }
}
