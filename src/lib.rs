use std::collections::{HashMap, HashSet};
use vector2math::*;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum Direction {
    UP,
    RIGHT,
    DOWN,
    LEFT,
}

/// A restricted direction is either UP or LEFT. Useful for NormalisedEdgePosition.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum RestrictedDirection {
    UP,
    LEFT,
}

const UP: (i32, i32) = (0, 1);
const RIGHT: (i32, i32) = (1, 0);
const DOWN: (i32, i32) = (0, -1);
const LEFT: (i32, i32) = (-1, 0);

type Glue = u8;

#[derive(Copy, Clone, Debug)]
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

type TilePosition = (i32, i32);
/// We denote edge by pairs such as ((0,0), LEFT)
/// The relation between TilePosition and EdgePosition is that
/// the right edge of the tile at TilePosition(0,0) is ((0,0), (0,-1)) = EdgePosition((0,0), LEFT).
type EdgePosition = ((i32, i32), Direction);
// Normalised edge position either has direction LEFT or UP.
type NormalisedEdgePosition = ((i32, i32), RestrictedDirection);

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

/// Returns the tile positions that are adjacent to the given edge position.
pub fn adjacent_tile_positions(edge_position: NormalisedEdgePosition) -> [TilePosition; 2] {
    match edge_position {
        (base_pos, RestrictedDirection::UP) => [base_pos, base_pos.add(RIGHT)],
        (base_pos, RestrictedDirection::LEFT) => [base_pos, base_pos.add(DOWN)],
    }
}

#[derive(Debug)]
pub enum TileTypeOrImpossible {
    TileType(TileType),
    NoTileTypeCanFit,
}

/// A tile assembly is a set of tiles positioned on a grid.
pub struct TileAssembly {
    tiles: HashMap<TilePosition, TileTypeOrImpossible>,
    edges: HashMap<NormalisedEdgePosition, Glue>,
    tileset: TileSet,
    /// The frontier is the set of tile positions that have at least one edge defined.
    current_frontier: HashSet<TilePosition>,
}

impl TileAssembly {
    pub fn solve(&mut self) {
        loop {
            let old_frontier = self.current_frontier.clone();
            self.solve_frontier();
            if old_frontier == self.current_frontier {
                break;
            }
        }
    }

    pub fn solve_frontier(&mut self) {
        let mut positions_to_remove: Vec<TilePosition> = Vec::new();
        for tile_position in self.current_frontier.clone().iter() {
            let adjacent_edges = self.adjacent_edges(tile_position);
            let matching_tiles = self.tileset.matching_tiles(&adjacent_edges);
            if matching_tiles.is_empty() {
                self.tiles
                    .insert(*tile_position, TileTypeOrImpossible::NoTileTypeCanFit);
            } else if matching_tiles.len() == 1 {
                let tile_type_index = matching_tiles[0];
                self.add_tile_from_tileset(tile_position, tile_type_index)
                    .unwrap();
                positions_to_remove.push(*tile_position);
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
                println!("Adding edge {:?}", edge_to_add);
                to_return.add_edge(&edge_to_add, 0).unwrap();
                current_position = current_position.add(LEFT);
            } else {
                let edge_to_add = (current_position, Direction::DOWN);
                println!("Adding edge {:?}", edge_to_add);
                to_return.add_edge(&edge_to_add, 1).unwrap();
                current_position = current_position.add(DOWN);

                let edge_to_add = (current_position, Direction::LEFT);
                println!("Adding edge {:?}", edge_to_add);
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
