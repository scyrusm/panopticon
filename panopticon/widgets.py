class Panopticopter:
    """ """

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
        """ """

        import numpy as np
        from bqplot import LinearScale, Scatter, Figure, Axis
        from ipywidgets import HBox, VBox, Button, widgets

        self.xs, self.ys = LinearScale(), LinearScale()
        self.button_sel0 = Button(description='Reset axes', icon='fa-refresh')
        self.button_sel1 = Button(description='Zoom', icon='fa-arrows')
        self.button_sel2 = Button(description='Select', icon='fa-crosshairs')
        self.button_sel3 = Button(description='Resubset', icon='fa-refresh')
        self.button_calc = Button(
            description='Calculate!',
            icon='fa-calculator',
        )
        self.button_layout = VBox([
            HBox([self.button_sel0, self.button_sel1, self.button_sel2]),
            HBox([self.button_sel3, self.button_calc])
        ])
        self.diffex_output = widgets.Output()
        self.scatter = Scatter(x=loom.ca[self.subset_mask][self.coord1],
                               y=loom.ca[self.subset_mask][self.coord2],
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
        """ """
        from bqplot.interacts import (PanZoom)
        panzoom = PanZoom(scales={"x": [self.xs], "y": [self.ys]})
        self.fig.interaction = panzoom

    def _reset(self):
        """ """
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
        """ """
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
            [True, False], p=[subset_p, 1 - subset_p], size=loom.shape[1])
        self.scatter.x = loom.ca[self.subset_mask][self.coord1]
        self.scatter.y = loom.ca[self.subset_mask][self.coord2]

    def _lasso_setup(self):
        """ """
        if self.fig.interaction is not None:
            self.fig.interaction.close()
            self.scatter.colors = ['gray']  #*self.subset_mask.sum()
            self.selection1 = []
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
        self.fig.interaction.color = 'black'

        # _which_to_select = self._which_to_select
        def make_selection():
            """ """
            self.scatter.unobserve_all()
            if 'blue' not in self.scatter.colors:
                which_to_select = 'blue'
            elif 'red' not in self.scatter.colors:
                which_to_select = 'red'
            else:
                which_to_select = None
            if which_to_select == "blue":

                self.selection1 = self.scatter.selected
                self._ignore_observe = True
                lasso_sel.reset()
                self._ignore_observe = False

                #recolor()
            elif which_to_select == "red":
                self.selection2 = self.scatter.selected
                self._ignore_observe = True
                lasso_sel.reset()
                self._ignore_observe = False

                #recolor()
            self.scatter.observe(select_and_recolor, names=['selected'])

        def select_and_recolor(self):
            """ """
            make_selection()
            recolor()

        def recolor():
            """ """
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

        #  lasso_sel.o
        def button0_click(widget):
            """

            Parameters
            ----------
            widget :
                

            Returns
            -------

            """
            self._reset()

        def button1_click(widget):
            """

            Parameters
            ----------
            widget :
                

            Returns
            -------

            """
            self._zoom_setup()

        def button2_click(widget):
            """

            Parameters
            ----------
            widget :
                

            Returns
            -------

            """
            self._lasso_setup()

        def button3_click(widget):
            """

            Parameters
            ----------
            widget :
                

            Returns
            -------

            """
            self.resubset()

        def button_calc_click(widget):
            """

            Parameters
            ----------
            widget :
                

            Returns
            -------

            """
            self.diffex_output.clear_output()
            from panopticon.analysis import get_cluster_differential_expression
            mask1 = np.isin(
                loom.ca['cellname'],
                loom.ca['cellname'][self.subset_mask][self.selection1])
            mask2 = np.isin(
                loom.ca['cellname'],
                loom.ca['cellname'][self.subset_mask][self.selection2])
            with self.diffex_output:

                if mask1 is not None and mask2 is not None:
                    diffex = get_cluster_differential_expression(
                        self.loom,
                        self.layername,
                        mask1=mask1,
                        mask2=mask2,
                        verbose=True)
                interesting_columns = [
                    'pvalue', 'CommonLanguageEffectSize', 'GeneAlternateName'
                ]
                diffex = diffex[interesting_columns]
                diffex.columns = ['p', 'CLES', 'gene']
                d1 = diffex.sort_values('CLES').head(10)
                d2 = diffex.sort_values('CLES', ascending=False).head(10)
                print('UPregulated in RED')
                display(diffex.sort_values('CLES').head(10), exclude='index')
                print('UPregulated in BLUE')
                display(diffex.sort_values('CLES', ascending=False).head(10))

        #  diffex_output1.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize').head(10).to_html()))
        #with diffex_output2:
        #  diffex_output2.append_display_data(HTML(diffex.sort_values('CommonLanguageEffectSize', ascending=False).head(10).to_html()))
        self.button_sel0.on_click(button0_click)
        self.button_sel1.on_click(button1_click, )
        self.button_sel2.on_click(button2_click)
        self.button_sel3.on_click(button3_click)
        self.button_calc.on_click(button_calc_click)

        from IPython.display import display
        #scatter_lasso.enable_move = True
        #fig_lasso
        final_layout = HBox(
            [VBox([self.fig, self.button_layout]), self.diffex_output])
        return final_layout
