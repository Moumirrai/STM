#import "@preview/lilaq:0.6.0" as lq

#let data = lq.load-txt(read("vornoi_fidelity.csv"))
#let (x_raw, E_mean, E_std, v_mean, v_std) = data

//#{ x_raw = x_raw.rev() }

// Data uz jsou serazena od hrube po jemnou, bez obraceni.
#let xs = x_raw

#let E_upper = E_mean.zip(E_std).map(((m, s)) => m + s)
#let E_lower = E_mean.zip(E_std).map(((m, s)) => m - s)
#let v_upper = v_mean.zip(v_std).map(((m, s)) => m + s)
#let v_lower = v_mean.zip(v_std).map(((m, s)) => m - s)

= Popis grafu
Předložené grafy zachycují vliv poloměru kontrolního obvodu používaného při
umisťování bodů, z nichž je následně generována Voronoiho teselační struktura.
Na takto vytvořené struktuře jsou určovány efektivní materiálové parametry,
konkrétně Youngův modul $E$ (horní graf) a Poissonovo číslo $nu$ (dolní graf).

Parametr poloměru kontrolního obvodu je měněn od hodnoty $0.6$
(hrubší struktura) po hodnotu $0.05$ (jemnější struktura).
Pro každou hodnotu je provedeno $300$ nezávislých iterací,
ze kterých je vyhodnocena střední hodnota a směrodatná odchylka.

Plná křivka v grafech reprezentuje střední hodnotu sledované veličiny,
zatímco barevné pásmo odpovídá intervalu střední hodnota ± směrodatná odchylka.
Výsledky ukazují, že střední hodnoty parametrů se v rámci zkoumaného rozsahu
mění pouze omezeně, zatímco se zjemňováním struktury se mění velikost rozptylu
výsledků.

Dolní ilustrace slouží k názornému porovnání hrubé a jemné varianty struktury.

#figure(
  grid(
    columns: 1,
    rows: 2,
    gutter: 1em,

    lq.diagram(
      xlabel: [Poloměr kontrolního obvodu],
      ylabel: [$E$ [Pa]],
      legend: (position: top + right),
      width: 10cm,
      height: 5cm,
      // Logaritmická osa X
      xscale: lq.scale.log(base: 10),
      xaxis: (inverted: true),
      lq.fill-between(
        xs,
        E_lower,
        y2: E_upper,
        fill: blue.lighten(60%),
        label: [$E$ směrodatná odchylka],
      ),
      lq.plot(
        xs,
        E_mean,
        mark: "o",
        mark-size: 3pt,
        color: blue,
        label: [$E$ střední hodnota],
      ),
    ),

    lq.diagram(
      xlabel: [Poloměr kontrolního obvodu],
      ylabel: [$nu$],
      legend: (position: bottom + right),
      xscale: lq.scale.log(base: 10),
      xaxis: (inverted: true),
      width: 10cm,
      height: 5cm,
      lq.fill-between(
        xs,
        v_lower,
        y2: v_upper,
        fill: red.lighten(60%),
        label: [$nu$ směrodatná odchylka],
      ),
      lq.plot(
        xs,
        v_mean,
        mark: "o",
        mark-size: 3pt,
        color: red,
        label: [$nu$ střední hodnota],
      ),
    ),
  ),
)
