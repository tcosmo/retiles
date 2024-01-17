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

pub struct TileSet {
    tiles: Vec<TileType>,
}

impl TileSet {
    /// Returns the Collatz tileset, cf. Tristan's thesis
    /// # Examples
    ///  
    /// ```
    /// use retiles::get_collatz_tileset;
    /// let tileset = get_collatz_tileset();
    /// ```
    pub fn get_collatz_tileset() -> Self {
        TileSet {
            tiles: vec![
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

/// The edge (0,0) -> (-1,0) can be defined either by ((0,0), LEFT) or ((-1,0), RIGHT).
/// The normalised edge position is the one where the direction is always LEFT or UP.
pub fn normalise_edge_position(edge_position: EdgePosition) -> NormalisedEdgePosition {
    let (base_pos, direction) = edge_position;
    match direction {
        Direction::UP => (base_pos, RestrictedDirection::UP),
        Direction::LEFT => (base_pos, RestrictedDirection::LEFT),
        Direction::RIGHT => (base_pos.add(RIGHT), RestrictedDirection::LEFT),
        Direction::DOWN => (base_pos.add(DOWN), RestrictedDirection::UP),
    }
}
pub struct AdjacentEdges {
    up: NormalisedEdgePosition,
    right: NormalisedEdgePosition,
    down: NormalisedEdgePosition,
    left: NormalisedEdgePosition,
}

/// Returns the tile positions that are adjacent to the given edge position.
pub fn adjacent_tile_positions(edge_position: NormalisedEdgePosition) -> [TilePosition; 2] {
    match edge_position {
        (base_pos, RestrictedDirection::UP) => [base_pos.add(LEFT), base_pos.add(RIGHT)],
        (base_pos, RestrictedDirection::LEFT) => [base_pos.add(UP), base_pos.add(DOWN)],
    }
}

/// Returns the edge positions that are adjacent to the given tile position in order UP, RIGHT, DOWN, LEFT.
pub fn adjacent_edge_positions(tile_position: TilePosition) -> AdjacentEdges {
    AdjacentEdges {
        up: (tile_position.add(UP), RestrictedDirection::LEFT),
        right: (tile_position, RestrictedDirection::UP),
        down: (tile_position, RestrictedDirection::LEFT),
        left: (tile_position.add(LEFT), RestrictedDirection::UP),
    }
}

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
    /// Adds an edge to the tile assembly.
    ///
    /// # Returns
    /// - Err(()) if the edge cannot be added (e.g. because it is already there)
    /// - OK(()) if the edge was added successfully
    pub fn add_edge(&mut self, edge_position: EdgePosition, glue: Glue) -> Result<(), ()> {
        let edge_position = normalise_edge_position(edge_position);
        match self.edges.get(&edge_position) {
            Some(_) => {
                return Result::Err(());
            }
            None => {
                self.edges.insert(edge_position, glue);
                for tile_position in adjacent_tile_positions(edge_position) {
                    match self.tiles.get(&tile_position) {
                        Some(TileTypeOrImpossible::TileType(_)) => {
                            // If there was a tile nearby the edge would have been added already.
                            assert!(false);
                        }
                        Some(TileTypeOrImpossible::NoTileTypeCanFit) => {
                            // This case where a side of the edge cannot fit any tile should be rare in practice but we allow it.
                            self.current_frontier.insert(tile_position);
                        }
                        None => {
                            self.current_frontier.insert(tile_position);
                        }
                    }
                }
                return Result::Ok(());
            }
        }
    }

    /// Returns true if the tile can be placed at the given position.
    /// This means that there is no tile at the position and that all
    /// potentially already defined edges match.
    ///
    pub fn is_tile_placeable_at_position(
        &self,
        tile_position: TilePosition,
        tile_type: TileType,
    ) -> bool {
        if self.tiles.contains_key(&tile_position) {
            return false;
        }
        let adjacent_edges = adjacent_edge_positions(tile_position);
        if let Some(&glue) = self.edges.get(&adjacent_edges.up) {
            if glue != tile_type.up {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.right) {
            if glue != tile_type.right {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.down) {
            if glue != tile_type.down {
                return false;
            }
        }
        if let Some(&glue) = self.edges.get(&adjacent_edges.left) {
            if glue != tile_type.left {
                return false;
            }
        }
        true
    }

    /// Add a tile from the tile set to the tile assembly.
    pub fn add_tile_from_tileset(
        &mut self,
        tile_position: TilePosition,
        tile_type_index: usize,
    ) -> Result<(), ()> {
        if tile_type_index >= self.tileset.tiles.len() {
            return Result::Err(());
        }
        let tile_type = self.tileset.tiles[tile_type_index];
        self.add_tile(tile_position, tile_type)
    }

    /// Add a tile to the tile  assembly.
    ///
    /// # Returns
    /// - Err(()) if the tile cannot be added (e.g. because it is already there or one of the edges that is already defined mismatches)
    /// - OK(()) if the tile was added successfully
    pub fn add_tile(&mut self, tile_position: TilePosition, tile_type: TileType) -> Result<(), ()> {
        if !self.is_tile_placeable_at_position(tile_position, tile_type) {
            return Result::Err(());
        }

        let adjacent_edges = adjacent_edge_positions(tile_position);

        self.edges.insert(adjacent_edges.up, tile_type.up);
        self.edges.insert(adjacent_edges.right, tile_type.right);
        self.edges.insert(adjacent_edges.down, tile_type.down);
        self.edges.insert(adjacent_edges.left, tile_type.left);

        self.tiles
            .insert(tile_position, TileTypeOrImpossible::TileType(tile_type));

        return Result::Ok(());
    }

    pub fn new(tileset: TileSet) -> Self {
        TileAssembly {
            tiles: HashMap::new(),
            edges: HashMap::new(),
            tileset: tileset,
            current_frontier: HashSet::new(),
        }
    }

    pub fn assembly_from_Collatz_parity_vector(parity_vector: Vec<u8>) -> Self {
        let mut to_return: TileAssembly = TileAssembly::new(TileSet::get_collatz_tileset());
        let mut current_position: TilePosition = (0, 0);
        for parity in parity_vector {
            if parity == 0 {
                to_return.add_edge((current_position,Direction::LEFT), 0).unwrap();
                current_position = current_position.add(LEFT);
            } else {
                to_return.add_edge((current_position,Direction::DOWN), 1).unwrap();
                current_position = current_position.add(LEFT);
                to_return.add_edge((current_position,Direction::LEFT), 0).unwrap();
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
        let ta = TileAssembly::assembly_from_Collatz_parity_vector(vec![0, 1, 0, 1, 1]);

        for (edge_position, glue) in ta.edges.iter() {
            println!("{:?} {:?}", edge_position, glue);
        }

        for position in ta.current_frontier {
            println!("{:?}", position);
        }

        assert_eq!(1, 1);
    }
}
