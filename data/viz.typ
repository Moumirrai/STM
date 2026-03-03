#import "voronoi.typ": points, connections, displacements
#import "@preview/cetz:0.4.2"

#let visualize_lite(
  structure,
  viz_scale: 1,
  displacement_scale: 1,
) = {
  cetz.canvas({
    import cetz.draw: *

    let scale = 1 / viz_scale

    set-viewport((0, 0), (1, 1), bounds: (scale, scale))

    let points = structure.points
    let base_scale_ratio = 0.14

    let displacedPoints = ()

    for (index, item) in structure.displacements.enumerate() {
      let original_point = points.at(index)
      let displaced_point = (
        original_point.at(0) + item.at(0) * displacement_scale * base_scale_ratio,
        original_point.at(1) + item.at(1) * displacement_scale * base_scale_ratio,
      )
      displacedPoints.push(displaced_point)
    }

    for item in structure.connections {
      on-layer(-1, {
        line(
          displacedPoints.at(int(item.at(0))),
          displacedPoints.at(int(item.at(1))),
          stroke: (paint: red),
        )
      })
      line(
        points.at(int(item.at(0))),
        points.at(int(item.at(1))),
        stroke: black,
      )
    }
  })
}

#align(center)[
  #visualize_lite(
    (
      points: points,
      connections: connections,
      displacements: displacements,
    ),
    viz_scale: 1,
    displacement_scale:3,
  )
]