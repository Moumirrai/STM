#import "@preview/lilaq:0.6.0" as lq

#let data = lq.load-txt(read("vornoi_fidelity.csv"))
#let (x_raw, E_mean, E_std, v_mean, v_std) = data

// Data jsou od hrubé (x=1) po jemnou (x=0.1) síť
// Použijeme kroky 1–20 jako osu X, takže jemnost roste doprava
#let steps = range(1, 21).map(i => float(i))

#let E_upper = E_mean.zip(E_std).map(((m, s)) => m + s)
#let E_lower = E_mean.zip(E_std).map(((m, s)) => m - s)
#let v_upper = v_mean.zip(v_std).map(((m, s)) => m + s)
#let v_lower = v_mean.zip(v_std).map(((m, s)) => m - s)

// Obrátíme pořadí dat (index 0 = nejhrubší, index 19 = nejjemnější)
#let rev(arr) = arr.rev()

#figure(
  grid(
    columns: 1,
    rows: 2,
    gutter: 1em,

    lq.diagram(
      xlabel: [Krok zjemnění (hrubá → jemná)],
      ylabel: [$E$ [Pa]],
      lq.fill-between(
        rev(steps),
        rev(E_lower),
        y2: rev(E_upper),
        fill: blue.lighten(60%),
        label: [$E$ pásmo $p m sigma$],
      ),
      lq.plot(
        rev(steps),
        rev(E_mean),
        mark: "o",
        mark-size: 3pt,
        color: blue,
        label: [$E$ střední hodnota],
      ),
    ),

    lq.diagram(
      xlabel: [Krok zjemnění (hrubá → jemná)],
      ylabel: [$nu$],
      lq.fill-between(
        rev(steps),
        rev(v_lower),
        y2: rev(v_upper),
        fill: red.lighten(60%),
        label: [$nu$ pásmo $p m sigma$],
      ),
      lq.plot(
        rev(steps),
        rev(v_mean),
        mark: "o",
        mark-size: 3pt,
        color: red,
        label: [$nu$ střední hodnota],
      ),
    ),
  )
)