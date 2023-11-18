import porepy as pp
import numpy as np


class RectangularDomainOrthogonalFractures2d(pp.ModelGeometry):
    def set_fracture_network(self) -> None:
        # Length scale:
        ls = 1 / self.units.m

        num_fracs = self.params.get("num_fracs", 1)
        domain = {"xmin": 0, "xmax": 2 * ls, "ymin": 0, "ymax": 1 * ls}
        if num_fracs == 0:
            p = np.zeros((2, 0), dtype=float) * ls
            e = np.zeros((2, 0), dtype=int)
        elif num_fracs == 1:
            p = np.array([[0, 2], [0.5, 0.5]]) * ls
            e = np.array([[0], [1]])
        elif num_fracs == 2:
            p = np.array([[0, 2, 0.5, 0.5], [1, 1, 0, 1]]) * ls
            e = np.array([[0, 2], [1, 3]])
        else:
            raise ValueError("Only 0, 1 or 2 fractures supported.")
        self.fracture_network = pp.FractureNetwork2d(p, e, domain)

    def mesh_arguments(self) -> dict:
        # Length scale:
        ls = 1 / self.units.m
        return {"mesh_size_frac": 0.2 * ls, "mesh_size_bound": 0.2 * ls}

class MassBalance(
    RectangularDomainOrthogonalFractures2d,
    pp.fluid_mass_balance.SinglePhaseFlow,
):
    ...

params = {}
#single_phase = pp.fluid_mass_balance.SinglePhaseFlow(params)
single_phase = MassBalance(params)
pp.run_time_dependent_model(single_phase, params)
pp.plot_grid(single_phase.mdg, "pressure")


class ModifiedFluidFlowBCs:
    def bc_values_darcy(self, subdomains: list[pp.Grid]) -> pp.ad.Array:
        """Dirichlet value p=1 on the west boundary of the domain.

        Parameters:
            subdomains: List of subdomains on which to define boundary conditions.

        Returns:
            Array of boundary values.

        """
        # Define boundary regions
        values = []
        for sd in subdomains:
            _, east, west, *_ = self.domain_boundary_sides(sd)
            val_loc = np.zeros(sd.num_faces)
            # See section on scaling for explanation of the conversion.
            val_loc[east] = self.fluid.convert_units(1, "Pa")
            val_loc[west] = self.fluid.convert_units(2, "Pa")
            values.append(val_loc)

        return pp.wrap_as_ad_array(np.hstack(values), name="bc_values_darcy")

    def bc_values_mobrho(self, subdomains: list[pp.Grid]) -> pp.ad.Array:
        """

        Parameters:
            subdomains: List of subdomains on which to define boundary conditions.

        Returns:
            Array of boundary values.

        """
        # Define boundary regions
        values = []
        for sd in subdomains:
            _, _, west, *_ = self.domain_boundary_sides(sd)
            val_loc = np.zeros(sd.num_faces)
            val_loc[west] = self.fluid.density() / self.fluid.viscosity()
            values.append(val_loc)
        return pp.wrap_as_ad_array(np.hstack(values), name="bc_values_mobrho")




class ModifiedSinglePhaseFlow(
    RectangularDomainOrthogonalFractures2d,
    ModifiedFluidFlowBCs,
    pp.fluid_mass_balance.SinglePhaseFlow,
):
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
            cell = sd.closest_cell(np.reshape([1.5, 0.8, 0], (3, 1)))
            val[cell] = 8
        return val


fluid_constants = pp.FluidConstants({"compressibility": 0.02, "viscosity": 0.1})
solid_constants = pp.SolidConstants({"permeability": 0.5, "porosity": 0.25})
material_constants = {"fluid": fluid_constants, "solid": solid_constants}
params = {"material_constants": material_constants}
single_phase = ModifiedSinglePhaseFlow(params)
pp.run_time_dependent_model(single_phase, params)
pp.plot_grid(single_phase.mdg, "pressure", figsize=(10, 8), linewidth=0, title="Pressure distribution")