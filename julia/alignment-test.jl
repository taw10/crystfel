using CrystFEL

panels = [Panel("q1", 512, 512, (-530, -530), 0.1, (1, 0, 0), (0, 1, 0), 100e-6, 1),
          Panel("q2", 512, 512, (28, -530),   0.1, (1, 0, 0), (0, 1, 0), 100e-6, 1),
          Panel("q3", 512, 512, (-530, 28),   0.1, (1, 0, 0), (0, 1, 0), 100e-6, 1),
          Panel("q4", 512, 512, (28, 28),     0.1, (1, 0, 0), (0, 1, 0), 100e-6, 1)]

detgeom = DetGeom(panels)

image = Image(detgeom)

cell = UnitCell(123, 45, 80, 90, 97, 90)

truth = predictreflections(image, cell)
