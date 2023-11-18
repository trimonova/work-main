import porepy as pp
import numpy as np

class ChangedGrid(pp.fluid_mass_balance.SinglePhaseFlow):
    """An Incompressible flow class with non-default mixed-dimensional grid.
    """

    def create_grid(self) -> None:
        """Create the grid bucket.

        A unit square grid with one vertical fracture extending from y=0.2 to
        y=0.5 is assigned.

        The method assigns the following attributes to self:
            mdg (pp.MixedDimensionalGrid): The produced mixed-dimensional grid.
            box (dict): The bounding box of the domain, defined through minimum and
                maximum values in each dimension.
        """
        # Use default mesh size if none are provided in the parameters passed to the class
        # on initialization
        mesh_args = self.params.get("mesh_args", {"mesh_size_frac": .1})
        endp = np.array([.2, .5])
        self.mdg, self.box = pp.md_grids_2d.single_vertical(mesh_args, endp, simplex=True)

params = {}
model_changed_grid = ChangedGrid(params)
pp.run_time_dependent_model(model_changed_grid, params)
pp.plot_grid(model_changed_grid.mdg, "pressure", figsize=[10,7])


class UnitaryFractureSource(ChangedGrid):
    """An IncompressibleFlow model with modified grid and
    unitary source term in the fracture."""

    def _source(self, sd: pp.Grid) -> np.ndarray:
        if sd.dim == self.mdg.dim_max():
            val = np.zeros(sd.num_cells)
        else:
            val = np.ones(sd.num_cells)
        return val


# We also prescribe a smaller mesh size:
params.update({"mesh_args": {"mesh_size_frac": 0.05}})
model_source = UnitaryFractureSource(params)
pp.run_time_dependent_model(model_source, params)
#pp.run_stationary_model(model_source, params)
pp.plot_grid(model_source.mdg, 'pressure', figsize=[10, 7])


class FractureAndMatrixSource(UnitaryFractureSource):
    """An IncompressibleFlow model with modified grid and
    source term in both matrix and fracture."""

    def _source(self, sd: pp.Grid) -> np.ndarray:
        """Return cell-wise source values.

        First call parent class method, then add value in matrix.

        Parameters:
            sd (Grid): The subdomain grid for which the source values are computed.

        Returns:
            np.ndarray (sd.num_cells): The cell-wise values.
        """
        # Gives 0 in matrix and 1 in fractures.
        val = super()._source(sd)
        # Add a source in the matrix cell at (0.7, 0.7):
        if sd.dim == self.mdg.dim_max():
            cell = sd.closest_cell(np.reshape([0.7, 0.7, 0], (3, 1)))
            val[cell] = 1
        return val


model = FractureAndMatrixSource(params)
#pp.run_stationary_model(model, params)
pp.run_time_dependent_model(model, params)
pp.plot_grid(model.mdg, 'pressure', figsize=[10, 7])


