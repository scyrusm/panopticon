class Panopticopter:

    def __init__(self,
                 loom,
                 layername,
                 coord1='log2(TP10k+1) PCA UMAP embedding 1',
                 coord2='log2(TP10k+1) PCA UMAP embedding 2',
                 n_points_to_display=1000):
        self.loom = loom
        self.layername = layername

        self.coord1 = coord1
        self.coord2 = coord2
        for coordi in [coord1, coord2]:
            if coordi not in loom.ca.keys():
                raise Exception("{} not in loom.ca.keys()".format(coordi))
        self.selection1 = []
        self.selection2 = []
        self.subset_mask = None
        self.fig = None
        self.n_point_to_display = n_points_to_display
        self.xs = None
        self.ys = None
        self.button_sel1 = None
        self.button_sel2 = None
        self.button_calc = None
        self.button_layout = None
        self.diffex_output = None
        self._which_to_select = "blue"
        self._ignore_observe = False
        import numpy as np
        subset_p = np.min([1, self.n_point_to_display / self.loom.shape[1]])
        self.subset_mask = np.random.choice([True, False],
                                            p=[subset_p, 1 - subset_p],
                                            size=loom.shape[1])

    def _fig_setup(self):

        import numpy as np
        from bqplot import LinearScale, Scatter, Figure, Axis
        from ipywidgets import HBox, VBox, Button, widgets

        self.xs, self.ys = LinearScale(), LinearScale()
        self.button_sel0 = Button(description='Reset axes', icon='fa-refresh')
        self.button_sel1 = Button(description='Zoom', icon='fa-arrows')
        self.button_sel2a = Button(description='Select', icon='fa-crosshairs')
        self.button_sel2a.style.button_color = 'blue'

        self.button_sel2b = Button(description='Select', icon='fa-crosshairs')
        self.button_sel2b.style.button_color = 'red'

        self.button_sel3 = Button(description='Resubset', icon='fa-refresh')
        self.button_calc = Button(
            description='Calculate!',
            icon='fa-calculator',
        )

        self.button_layout = VBox([
            HBox([self.button_sel0, self.button_sel1, self.button_sel2a]),
            HBox([self.button_sel3, self.button_calc, self.button_sel2b])
        ])
        self.diffex_output = widgets.Output()
        self.scatter = Scatter(x=self.loom.ca[self.subset_mask][self.coord1],
                               y=self.loom.ca[self.subset_mask][self.coord2],
                               scales={
                                   "x": self.xs,
                                   "y": self.ys
                               },
                               default_size=1,
                               colors=['gray'])
        xax, yax = Axis(scale=self.xs,
                        label="UMAP 1"), Axis(scale=self.ys,
                                              label="UMAP 2",
                                              orientation="vertical")
        self.fig = Figure(
            marks=[
                self.scatter,
            ],
            axes=[xax, yax],
            title="Panopticopter",
            #interaction=lasso_sel,
        )

    def _zoom_setup(self):
        from bqplot.interacts import (PanZoom)
        panzoom = PanZoom(scales={"x": [self.xs], "y": [self.ys]})
        self.fig.interaction = panzoom

    def _reset(self):
        import numpy as np
        xmin = float(np.min(self.loom.ca[self.coord1]))
        xmax = float(np.max(self.loom.ca[self.coord2]))
        ymin = float(np.min(self.loom.ca[self.coord1]))
        ymax = float(np.max(self.loom.ca[self.coord2]))
        self.scatter.scales['x'].min = xmin
        self.scatter.scales['x'].max = xmax
        self.scatter.scales['y'].min = ymin
        self.scatter.scales['y'].max = ymax

    def resubset(self):
        import numpy as np
        mask = np.array([True] * self.loom.shape[1])
        if self.xs.min is not None:
            mask *= self.loom.ca[self.coord1] >= self.xs.min
        if self.xs.max is not None:
            mask *= self.loom.ca[self.coord1] <= self.xs.max
        if self.ys.min is not None:
            mask *= self.loom.ca[self.coord2] >= self.ys.min
        if self.ys.max is not None:
            mask *= self.loom.ca[self.coord2] <= self.ys.max
        subset_p = np.min([1, self.n_point_to_display / mask.sum()])
        self.subset_mask = mask * np.random.choice(
            [True, False], p=[subset_p, 1 - subset_p], size=self.loom.shape[1])
        self.scatter.x = self.loom.ca[self.subset_mask][self.coord1]
        self.scatter.y = self.loom.ca[self.subset_mask][self.coord2]

    def _lasso_setup(self, which_to_select):
        if self.fig.interaction is not None:
            self.fig.interaction.close()
            if which_to_select == 'blue':
                self.scatter.colors = [
                    x.replace('blue', 'gray') for x in self.scatter.colors
                ]  #*self.subset_mask.sum()
                self.selection1 = []
            elif which_to_select == 'red':
                self.scatter.colors = [
                    x.replace('red', 'gray') for x in self.scatter.colors
                ]  #*self.subset_mask.sum()
                self.selection2 = []

        from bqplot.interacts import (
            LassoSelector,
            #PanZoom,
        )
        lasso_sel = LassoSelector()
        lasso_sel.marks = [
            self.scatter,
        ]
        self.fig.interaction = lasso_sel
        self.fig.interaction.color = which_to_select

        # _which_to_select = self._which_to_select
        def make_selection():
            self.scatter.unobserve_all()
            #if 'blue' not in self.scatter.colors:
            #    which_to_select = 'blue'
            #elif 'red' not in self.scatter.colors:
            #    which_to_select = 'red'
            #else:
            #    which_to_select = None
            if which_to_select == "blue":
                self.selection1 = self.scatter.selected
                #lasso_sel.reset()

                #recolor()
            elif which_to_select == "red":
                self.selection2 = self.scatter.selected
                #lasso_sel.reset()

                #recolor()
            self.scatter.observe(select_and_recolor, names=['selected'])

        def select_and_recolor(self):
            make_selection()
            recolor()

        def recolor():
            self.scatter.unobserve_all()
            if self.selection1 is None:
                self.selection1 = []
            if self.selection2 is None:
                self.selection2 = []
            colors = []
            for i in range(self.subset_mask.sum()):
                if i in self.selection1:
                    colors.append('blue')
                elif i in self.selection2:
                    colors.append('red')
                else:
                    colors.append('gray')
            self.scatter.colors = colors
            self.scatter.observe(select_and_recolor, names=['selected'])

        self.scatter.unobserve_all()
        self.scatter.observe(select_and_recolor, names=['selected'])

    def __call__(self):
        import numpy as np
        from bqplot import pyplot as plt
        #from bqplot import DateScale, LinearScale, Axis, Lines, Scatter, Bars, Hist, Figure
        from bqplot import LinearScale, Scatter, Figure, Axis
        from ipywidgets import HBox, VBox, Button, widgets

        self._fig_setup()
        self.selection1 = []
        self.selection2 = []

        #  lasso_sel.o
        def button0_click(widget):
            self._reset()

        def button1_click(widget):
            self._zoom_setup()

    # def button2_click(widget):
    #    self._lasso_setup()

        def button2a_click(widget):
            self._lasso_setup('blue')

        def button2b_click(widget):
            self._lasso_setup('red')

        def button3_click(widget):
            self.resubset()

        def button_calc_click(widget):
            self.diffex_output.clear_output()
            from panopticon.analysis import get_cluster_differential_expression
            import pandas as pd
            from IPython.display import HTML

            mask1 = np.isin(
                self.loom.ca['cellname'],
                self.loom.ca['cellname'][self.subset_mask][self.selection1])
            mask2 = np.isin(
                self.loom.ca['cellname'],
                self.loom.ca['cellname'][self.subset_mask][self.selection2])
            with self.diffex_output:

                if mask1 is not None and mask2 is not None:
                    diffex = get_cluster_differential_expression(
                        self.loom,
                        self.layername,
                        mask1=mask1,
                        mask2=mask2,
                        verbose=True)
                self.diffex = diffex
                if 'GeneAlternateName' in diffex.columns:
                    interesting_columns = [
                        'pvalue', 'CommonLanguageEffectSize',
                        'GeneAlternateName'
                    ]
                else:
                    interesting_columns = [
                        'pvalue', 'CommonLanguageEffectSize', 'gene'
                    ]
                diffex = diffex[interesting_columns]
                diffex.columns = ['p', 'CLES', 'gene']
                d1 = diffex.sort_values(
                    'CLES', ascending=False).head(10).reset_index(drop=True)
                d2 = diffex.sort_values(
                    'CLES', ascending=True).head(10).reset_index(drop=True)
                d1.columns = [x + '_blue' for x in d1.columns]
                d2.columns = [x + '_red' for x in d2.columns]

                d1d2 = pd.concat([d1, d2], axis=1)
                self.d1d2 = d1d2

                display(
                    d1d2.style.set_properties(**{
                        'background-color': 'blue'
                    },
                                              subset=d1.columns).
                    set_properties(**{
                        'background-color': 'red'
                    },
                                   subset=d2.columns).set_properties(
                                       **{'color': 'white'}))
            # print(HTML(d1d2.to_html(index=False)))

        #  diffex_output1.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize').head(10).to_html()))
        #with diffex_output2:
        #  diffex_output2.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize', ascending=False).head(10).to_html()))
        self.button_sel0.on_click(button0_click)
        self.button_sel1.on_click(button1_click, )
        self.button_sel2a.on_click(button2a_click)
        self.button_sel2b.on_click(button2b_click)

        self.button_sel3.on_click(button3_click)
        self.button_calc.on_click(button_calc_click)

        from IPython.display import display
        #scatter_lasso.enable_move = True
        #fig_lasso
        final_layout = HBox(
            [VBox([self.fig, self.button_layout]), self.diffex_output])
        return final_layout
