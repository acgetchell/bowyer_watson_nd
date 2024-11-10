use std::collections::HashSet;
use ordered_float::OrderedFloat;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct Point {
    coords: Vec<OrderedFloat<f64>>,
}

impl Point {
    fn dimension(&self) -> usize {
        self.coords.len()
    }
}

#[derive(Debug, Clone, PartialEq)]
struct Simplex {
    vertices: Vec<Point>,
}

impl Simplex {
    /// Checks if the point is inside the circumsphere of the simplex.
    fn contains_in_circumsphere(&self, point: &Point) -> bool {
        let n = self.vertices.len() - 1; // Dimension of the simplex

        // Build the augmented matrix for the determinant calculation
        let mut matrix = vec![vec![0.0; n + 2]; n + 2];

        for i in 0..=n {
            for j in 0..n {
                matrix[i][j] = (self.vertices[i].coords[j] - point.coords[j]).into();
            }
            matrix[i][n] = self.vertices[i]
                .coords
                .iter()
                .zip(&point.coords)
                .map(|(vi, pi)| (vi - pi).powi(2))
                .sum();
            matrix[i][n + 1] = 1.0;
        }

        // Last row
        for j in 0..=n + 1 {
            matrix[n + 1][j] = if j == n { 0.0 } else { 1.0 };
        }

        // Calculate the determinant
        let det = determinant(&matrix);

        // For the Delaunay triangulation, we check if the determinant is positive
        det > 0.0
    }

    /// Checks if the simplex contains the given point as a vertex.
    fn contains_vertex(&self, point: &Point) -> bool {
        self.vertices.contains(point)
    }

    /// Returns all faces (subsimplices of dimension n-1) of the simplex.
    fn faces(&self) -> Vec<Simplex> {
        let mut faces = Vec::new();
        let n = self.vertices.len();
        for i in 0..n {
            let mut face_vertices = self.vertices.clone();
            face_vertices.remove(i);
            faces.push(Simplex {
                vertices: face_vertices,
            });
        }
        faces
    }
}

/// Calculates the determinant of a square matrix.
fn determinant(matrix: &[Vec<f64>]) -> f64 {
    let n = matrix.len();
    if n == 1 {
        return matrix[0][0];
    }
    let mut det = 0.0;
    for i in 0..n {
        let mut sub_matrix = Vec::new();
        for j in 1..n {
            let mut row = Vec::new();
            for k in 0..n {
                if k != i {
                    row.push(matrix[j][k]);
                }
            }
            sub_matrix.push(row);
        }
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        det += sign * matrix[0][i] * determinant(&sub_matrix);
    }
    det
}

fn bowyer_watson(points: &[Point]) -> Vec<Simplex> {
    let dimension = points[0].dimension();

    // Create a super-simplex that encompasses all the points
    let mut max_coords = vec![f64::NEG_INFINITY; dimension];
    let mut min_coords = vec![f64::INFINITY; dimension];

    // Find the bounding box of all points
    for p in points {
        for i in 0..dimension {
            if p.coords[i] < ordered_float::OrderedFloat(min_coords[i]) {
                min_coords[i] = *p.coords[i];
            }
            if p.coords[i] > ordered_float::OrderedFloat(max_coords[i]) {
                max_coords[i] = *p.coords[i];
            }
        }
    }

    let mut delta = vec![0.0; dimension];
    for i in 0..dimension {
        delta[i] = max_coords[i] - min_coords[i];
        if delta[i] == 0.0 {
            delta[i] = 1.0; // Avoid zero delta
        }
    }

    // Construct super-simplex vertices
    let mut super_vertices = Vec::new();
    for i in 0..=dimension {
        let mut coords = Vec::new();
        for j in 0..dimension {
            let mut value = min_coords[j] - 10.0 * delta[j];
            if i == j {
                value = max_coords[j] + 10.0 * delta[j];
            }
            coords.push(OrderedFloat(value));
        }
        super_vertices.push(Point { coords });
    }

    // Initialize the triangulation with the super-simplex
    let super_simplex = Simplex {
        vertices: super_vertices.clone(),
    };
    let mut triangulation = Vec::new();
    triangulation.push(super_simplex);

    // Iterate over each point to construct the triangulation
    for p in points {
        let mut bad_simplices = Vec::new();

        // Find all simplices that are invalidated by the point insertion
        for s in &triangulation {
            if s.contains_in_circumsphere(p) {
                bad_simplices.push(s.clone());
            }
        }

        // Find the boundary (polytope hole) created by the bad simplices
        let mut boundary_faces = Vec::new();
        for bs in &bad_simplices {
            for face in bs.faces() {
                let mut is_shared = false;
                for bs2 in &bad_simplices {
                    if bs == bs2 {
                        continue;
                    }
                    if bs2.faces().contains(&face) {
                        is_shared = true;
                        break;
                    }
                }
                if !is_shared {
                    boundary_faces.push(face);
                }
            }
        }

        // Remove the bad simplices from the triangulation
        triangulation.retain(|s| !bad_simplices.contains(s));

        // Re-triangulate the polytope hole with new simplices formed by the point
        for face in boundary_faces {
            let mut new_simplex_vertices = face.vertices.clone();
            new_simplex_vertices.push(p.clone());
            triangulation.push(Simplex {
                vertices: new_simplex_vertices,
            });
        }
    }

    // Remove simplices that include the super-simplex's vertices
    let super_vertices_set: HashSet<Point> = super_vertices.into_iter().collect();
    triangulation.retain(|s| {
        !s.vertices
            .iter()
            .any(|v| super_vertices_set.contains(v))
    });

    triangulation
}

fn main() {
    // Example in 3D space
    let points = vec![
        Point {
            coords: vec![ordered_float::OrderedFloat(0.2), ordered_float::OrderedFloat(0.4), ordered_float::OrderedFloat(0.1)],
        },
        Point {
            coords: vec![ordered_float::OrderedFloat(0.5), ordered_float::OrderedFloat(0.1), ordered_float::OrderedFloat(0.3)],
        },
        Point {
            coords: vec![ordered_float::OrderedFloat(0.9), ordered_float::OrderedFloat(0.7), ordered_float::OrderedFloat(0.2)],
        },
        Point {
            coords: vec![ordered_float::OrderedFloat(0.4), ordered_float::OrderedFloat(0.9), ordered_float::OrderedFloat(0.5)],
        },
        Point {
            coords: vec![ordered_float::OrderedFloat(0.7), ordered_float::OrderedFloat(0.3), ordered_float::OrderedFloat(0.8)],
        },
    ];

    let simplices = bowyer_watson(&points);

    // Print the resulting triangulation
    for (i, s) in simplices.iter().enumerate() {
        println!("Simplex {}:", i + 1);
        for v in &s.vertices {
            println!("  Point: {:?}", v.coords);
        }
        println!();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determinant() {
        let matrix = vec![
            vec![1.0, 2.0, 3.0],
            vec![4.0, 5.0, 6.0],
            vec![7.0, 8.0, 9.0],
        ];
        assert_eq!(determinant(&matrix), 0.0);

        let matrix = vec![
            vec![1.0, 2.0],
            vec![3.0, 4.0],
        ];
        assert_eq!(determinant(&matrix), -2.0);

        let matrix = vec![
            vec![2.0],
        ];
        assert_eq!(determinant(&matrix), 2.0);
    }

    #[test]
    fn test_contains_vertex() {
        let point = Point { coords: vec![ordered_float::OrderedFloat(1.0), ordered_float::OrderedFloat(2.0)] };
        let simplex = Simplex { vertices: vec![point.clone()] };
        assert!(simplex.contains_vertex(&point));

        let other_point = Point { coords: vec![ordered_float::OrderedFloat(3.0), ordered_float::OrderedFloat(4.0)] };
        assert!(!simplex.contains_vertex(&other_point));
    }

    #[test]
    fn test_faces() {
        let point1 = Point { coords: vec![ordered_float::OrderedFloat(1.0), ordered_float::OrderedFloat(2.0)] };
        let point2 = Point { coords: vec![ordered_float::OrderedFloat(3.0), ordered_float::OrderedFloat(4.0)] };
        let point3 = Point { coords: vec![ordered_float::OrderedFloat(5.0), ordered_float::OrderedFloat(6.0)] };
        let simplex = Simplex { vertices: vec![point1.clone(), point2.clone(), point3.clone()] };

        let faces = simplex.faces();
        assert_eq!(faces.len(), 3);
        assert!(faces.contains(&Simplex { vertices: vec![point2.clone(), point3.clone()] }));
        assert!(faces.contains(&Simplex { vertices: vec![point1.clone(), point3.clone()] }));
        assert!(faces.contains(&Simplex { vertices: vec![point1.clone(), point2.clone()] }));
    }

    #[test]
    fn test_bowyer_watson() {
        let points = vec![
            Point { coords: vec![ordered_float::OrderedFloat(0.0), ordered_float::OrderedFloat(0.0)] },
            Point { coords: vec![ordered_float::OrderedFloat(1.0), ordered_float::OrderedFloat(0.0)] },
            Point { coords: vec![ordered_float::OrderedFloat(0.0), ordered_float::OrderedFloat(1.0)] },
            Point { coords: vec![ordered_float::OrderedFloat(1.0), ordered_float::OrderedFloat(1.0)] },
        ];

        let simplices = bowyer_watson(&points);
        assert!(!simplices.is_empty());
    }
}