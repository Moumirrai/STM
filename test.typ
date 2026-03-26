#import "@preview/lilaq:0.5.0" as lq

#let (x, E_mean, E_var, v_mean, v_var) = lq.load-txt(read("vornoi_fidelity.csv"))

//#let x = lq.linspace(0.001, 0.999, num: 100)

#figure(
  caption: "1x1 Bowtie Modul pružnosti",
)[
  #lq.diagram(
    legend: (position: (100% + .5em, 0%)),
    lq.plot(x, E_mean, mark: none, smooth: false, label: $E_"mean"$),
    lq.plot(x, E_mean.map(x => x * -1), mark: none, smooth: false, label: $E_"mean"$),
    xlabel: [Úhel $alpha$ (°)], ylabel: "Modul pružnosti (GPa)"
  )
]

#figure(
  caption: "1x1 Bowtie Modul pružnosti",
)[
  #lq.diagram(
    legend: (position: (100% + .5em, 0%)),
    lq.plot(x, E_var, mark: none, smooth: false, label: $E_"var"$),
    lq.plot(x, E_var.map(x => x * -1), mark: none, smooth: false, label: $E_"var"$),
    xlabel: [Úhel $alpha$ (°)], ylabel: "Modul pružnosti (GPa)"
  )
]