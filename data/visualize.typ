#import "@preview/cetz:0.4.1": canvas, draw, tree

#set page(
  paper: "a4",
  margin: 1.6cm,
)


= Debug

// set page format to a4 and margin to 1cm


#let data = json("line.json")

#canvas(length: 2cm, {
  import draw: *
  set-style(
    mark: (fill: black, scale: 2),
    stroke: (thickness: 5pt, cap: "round"),
    angle: (
      radius: 0.1,
      label-radius: .1,
      fill: green.lighten(80%),
      stroke: (paint: green.darken(50%)),
    ),
    content: (padding: 1pt),
  )

  let point(coords) = circle(coords, radius: 0.01, fill: black, stroke: none) //styl pro uzel

  for i in range(data.nodes.len()) { //pro každý uzel
    let node = data.nodes.at(i)
    point((float(eval(node.dx)), float(eval(node.dy)) * -1)) //tečka
    anchor(str(i + 1), (float(eval(node.dx)), float(eval(node.dy)) * -1)) //kotva pro text
    content(str(i + 1), anchor: "north", text(black)[$ #(i) $], padding: .1) //text (číslo uzlu)
  }

  for i in range(data.elements.len()) { //pro každý prut
    let element = data.elements.at(i)
    let startNode = data.nodes.at(int(element.starting_node))
    let endNode = data.nodes.at(int(element.ending_node))
    let startCoords = (float(eval(startNode.dx)), float(eval(startNode.dy)) * -1)
    let endCoords = (float(eval(endNode.dx)), float(eval(endNode.dy)) * -1)
    line(startCoords, endCoords, stroke: black) //čára
    anchor("el" + str(i + 1), ((startCoords.at(0) + endCoords.at(0)) / 2, (startCoords.at(1) + endCoords.at(1)) / 2)) //kotva do středu prutu
    content("el" + str(i + 1), anchor: "south-east", text(green)[$ #(i + 1) $], padding: .1) //text (číslo prutu)
  }
})

Total constraints: 4

Free DOFs: 6

Dependent DOFs: 0

Fixed DOFs: 4

#let ordered_dofs = (2, 3, 4, 5, 6, 7, 0, 1, 8, 9)

#import "@preview/teig:0.1.0": eigenvalues

#let kmat = (
  (1503600, -604800, -1050000, 0, -453600, 604800, 0, 0, 0, 0),
  (-604800, 806400, 0, 0, 604800, -806400, 0, 0, 0, 0),
  (-1050000, 0, 1503600, 604800, -453600, -604800, 0, 0, 0, 0),
  (0, 0, 604800, 2381400, -604800, -806400, 0, -1575000, 0, 0),
  (-453600, 604800, -453600, -604800, 3460800, -604800, -2100000, 0, -453600, 604800),
  (604800, -806400, -604800, -806400, -604800, 2419200, 0, 0, 604800, -806400),
  (0, 0, 0, 0, -2100000, 0, 2100000, 0, 0, 0),
  (0, 0, 0, -1575000, 0, 0, 0, 3150000, 0, -1575000),
  (0, 0, 0, 0, -453600, 604800, 0, 0, 453600, -604800),
  (0, 0, 0, 0, 604800, -806400, 0, -1575000, -604800, 2381400),
)

#import "@preview/pavemat:0.2.0": pavemat

#set math.mat(row-gap: 0.25em, column-gap: 0.1em)
#set text(size: 1em)

#let getDofName = dof => {
  let nodeid = calc.quo(dof, 2) + 1  // Integer division
  let direction = calc.rem(dof, 2)
  if (direction == 0) {
    return [$u_#nodeid$]
  } else {
    return [$w_#nodeid$]
  }
}

== Matrix K
#table(
  columns: 2,
  rows: 2,
  stroke: none,
  [],
  [
    #for i in range(ordered_dofs.len()) {
      if (i == 0) [#h(1.9em)#getDofName(ordered_dofs.at(i))] else [#h(3.08em)#getDofName(ordered_dofs.at(i))]
    }

  ],

  [
    //change line spacing
    #set par(leading: 0.57em)
    #for i in range(ordered_dofs.len()) {
      if (i == 0) [#v(0.8mm)#getDofName(ordered_dofs.at(i))\ ] else [#getDofName(ordered_dofs.at(i))\ ]
    }
  ],
  [
    #pad(left: -3mm)[
      #pavemat(
      pave: "SSSSSSDDDDDDWWWWWWAAAAAAdddddd(paint: red, thickness: 1pt)DDDDSSSSSSAAAAWWWWWW",
      stroke: (dash: "dashed", thickness: 1pt, paint: blue),
      fills: (
        "1-1": blue.transparentize(80%),
        "1-6": red.transparentize(80%),
      ),
      delim: "[",
    )[$mat(..kmat)$]
    ]
  ],
)

== Loads vector
\ 
#text(size: 1.4em)[
  $f_"free" = vec(u_2, w_2, u_3, w_3, u_4, w_4, delim: #none)vec(0, 0, 9, 4, 0, 0, delim: "[")$
]

== R vectors -> zero
