#! /usr/bin/env python

from pyschism.outputs import (
    PlotOutputCombined,
    PlotOutputStations,
    StationsOutput,
)


class OutputPlotCli:

    def __init__(self, args):

        if not args.resource.exists():
            raise OSError(f'ERROR: Invalid path {str(args.resource)}. '
                          'It does not exist.')

        self._args = args

        # args.output_type is at current a required argument but in the future
        # it doesn't have to be required.
        if args.output_type is None:
            # will never hit this part unless we remove the output_type req.
            raise NotImplementedError('Should plot everything.')

        else:
            getattr(self, args.output_type)()

    def surface(self):

        if self._args.resource.is_dir():
            raise NotImplementedError(
                'Uncombined outputs not yet supported, please pass a combined '
                'output and stay tuned for updates. This will be supported in '
                'the future.')

        elif self._args.resource.is_file():
            plotter = PlotOutputCombined(self._args.resource)
            plotter.plot(self._args.variable, show=True)

    def stations(self):
        plotter = PlotOutputStations(StationsOutput(self._args.resource))
        plotter.plot(self._args.variable)
