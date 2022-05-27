def panopticopter(loom,
        layername,
        coord=None,
        coord1=None,
        coord2=None):
    import numpy as np
    from panopticon.utilities import import_check
    exit_code = import_check("bqplot", 'pip install bqplot')
    if exit_code != 0:
        return
    exit_code = import_check("ipywidgets", 'pip install ipywidgets')
    if exit_code != 0:
        return
    if coord is not None and coord1 is None and coord2 is None:
        coord1 = coord+' 1'
        coord2 = coord+' 2'
    elif coord is None and coord1 is not None and coord2 is not None:
        pass
    else:
        raise Exception("Either coord must be None, xor coord1 and coord2 must be None")
    for coordi in [coord1, coord2]:
        if coordi not in loom.ca.keys():
            raise Exception("{} not in loom.ca.keys()".format(coordi))

    from bqplot import pyplot as plt

    from bqplot import DateScale, LinearScale, Axis, Lines, Scatter, Bars, Hist, Figure
    from bqplot.interacts import (
        FastIntervalSelector,
        IndexSelector,
        BrushIntervalSelector,
        BrushSelector,
        MultiSelector,
        LassoSelector,
        PanZoom,
        HandDraw,
    )
    from traitlets import link

    from ipywidgets import ToggleButtons, VBox, HTML, widgets

    lasso_sel = LassoSelector()
    subset_mask = np.random.choice([True, False],
                                   p=[0.05, 0.95],
                                   size=loom.shape[1])
    xs, ys = LinearScale(), LinearScale()
    #data = np.arange(20)
    from ipywidgets import FloatSlider, HBox, VBox, Button

    button_sel1 = Button(description='Selection 1')
    button_sel2 = Button(description='Selection 2')
    button_calc = Button(description='Calculate!')
    button_layout = HBox([button_sel1, button_sel2, button_calc])
    diffex_output = widgets.Output()
    #diffex_output2 = widgets.Output()

    scatter_lasso = Scatter(
        x=loom.ca[coord1][subset_mask],
        y=loom.ca[coord2][subset_mask],
        scales={
            "x": xs,
            "y": ys
        },
        default_size=1,
        colors=['gray'])
    #scatter_lasso = Scatter(x=data, y=data, scales={"x": xs, "y": ys}, colors=["skyblue"])
    #bar_lasso = Bars(x=data, y=data / 2.0, scales={"x": xs, "y": ys})
    xax_lasso, yax_lasso = Axis(scale=xs,
                                label="UMAP 1"), Axis(scale=ys,
                                                      label="UMAP 2",
                                                      orientation="vertical")
    fig_lasso = Figure(
        marks=[
            scatter_lasso,
        ],
        axes=[xax_lasso, yax_lasso],
        title="Panopticopter",
        interaction=lasso_sel,
    )
    lasso_sel.marks = [
        scatter_lasso,
    ]
    #fig_lasso.marks.set_metadata(1)
    selection1 = None
    selection2 = None

    def recolor(selection1, selection2):
        colors = []
        if selection1 is None:
            selection1 = []
        if selection2 is None:
            selection2 = []
        for i in range(subset_mask.sum()):
            if i in selection1:
                colors.append('blue')
            elif i in selection2:
                colors.append('red')
            else:
                colors.append('gray')
        #global scatter_lasso

        scatter_lasso.colors = colors
        #lasso_sel.marks = [scatter_lasso, ]

    def button1_click(widget):
        global selection1
        global selection2
        selection1 = scatter_lasso.selected
        lasso_sel.reset()
        recolor(selection1, selection2)

    def button2_click(widget):
        global selection1
        global selection2
        selection2 = scatter_lasso.selected
        lasso_sel.reset()
        recolor(selection1, selection2)

    def button_calc_click(widget):
        diffex_output.clear_output()
        from panopticon.analysis import get_cluster_differential_expression
        mask1 = np.isin(loom.ca['cellname'],
                        loom.ca['cellname'][subset_mask][selection1])
        mask2 = np.isin(loom.ca['cellname'],
                        loom.ca['cellname'][subset_mask][selection2])
        with diffex_output:

            if mask1 is not None and mask2 is not None:
                diffex = get_cluster_differential_expression(loom,
                                                            layername,
                                                             mask1=mask1,
                                                             mask2=mask2,
                                                             verbose=True)
            interesting_columns = [
                'pvalue', 'CommonLanguageEffectSize', 'GeneAlternateName'
            ]
            diffex = diffex[interesting_columns]
            diffex.columns = ['p', 'CLES', 'gene']
            print('UPregulated in RED')
            display(diffex.sort_values('CLES').head(10))
            print('UPregulated in BLUE')
            display(diffex.sort_values('CLES', ascending=False).head(10))

    #  diffex_output1.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize').head(10).to_html()))
    #with diffex_output2:
    #  diffex_output2.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize', ascending=False).head(10).to_html()))

    button_sel1.on_click(button1_click, )
    button_sel2.on_click(button2_click)
    button_calc.on_click(button_calc_click)

    from IPython.display import display
    #scatter_lasso.enable_move = True
    #fig_lasso
    final_layout = HBox([VBox([fig_lasso, button_layout]), diffex_output])
    return final_layout
