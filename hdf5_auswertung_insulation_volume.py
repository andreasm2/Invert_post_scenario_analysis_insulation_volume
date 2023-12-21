# -*- coding: iso-8859-15 -*-
"""

@author: Andreas Müller


Postanalysis module started by the Invert model
The aim of this module is to derive the calculated insulation volume and
windows exchanged for a given scenario in order to provide an estimate of material use for
a follow-up embedded energy and CO2-emission calculations
"""
import os
import sys
import time

import h5py

import numpy as np


def CreateInsulationData(SCENARIO_INFO, path_):
    """

    Calculates Insulation Demand per Simulation Period"""

    print("\n\n#####################")
    print("#####################")
    print("Calculate Insulation Demand")
    print("\n\n#####################")
    print("#####################\n\n")

    #########################
    #
    #
    # Some Calculation Assumptions
    #
    #
    #########################

    lambda_insulation = 0.04

    # Lower bound for calculated insulation thickness that needs to be exceeded in order to increase
    # thickness to min_ins_thickness_wall/min_ins_thickness_ceil/min_ins_thickness_floor
    min_calc_ins_thickness_to_increase = 0.03  # Meters

    # if min_calc_ins_thickness_to_increase is exceeded then
    # thickness must have at least specific thickness
    min_ins_thickness_wall = 0.08  # Meters
    min_ins_thickness_ceil = 0.08  # Meters
    min_ins_thickness_floor = 0.08  # Meters

    # upper boundary for insulation thickness
    max_ins_thickness_wall = 0.3  # Meters
    max_ins_thickness_ceil = 0.35  # Meters
    max_ins_thickness_floor = 0.3  # Meters

    # lower boundary for U-Value
    min_U_value = 0.08  # W/m2K

    # additional wall area resulting from buildings not being a QUADER
    vertical_and_horizontal_building_ABGESETZTHEIT = 1.2

    # additional factor (all Surfaces) resulting from Non-Optimimal distributed insulation
    # (e.g. Dachgiebel, Waermebruecken, Ecken, etc.
    # => in order to reach a given U-value, the average insulation thickness must be higher
    non_optimality_factor = 1.15

    # Estimated Demand for internal walls and floors
    # Share of buildings that need insulation material for internal walls and floors
    share_retrofit_floors = 0.45
    share_retrofit_walls = 0.25
    share_new_buildings_floors = 0.75
    share_new_buildings_walls = 0.5

    # internal Walls -> Assume Simple Cross => length + width of building
    thickness_internal_wall = 0.08
    thickness_internal_floor = 0.05

    #########################
    #
    #
    #
    #########################

    outputpath_split = path_.replace("\\", "/").replace("//", "/").split("/")
    dir_ = ("/").join(outputpath_split[:-1]) + "/"
    runnumberstring = outputpath_split[-1]

    output_path_ = dir_ + "ADD_RESULTS/INSULATION_DEMAND/"
    scenario_name = (path_.split("/")[-2]).replace("_", "\\_")

    filename = "%s/%s%s.hdf5" % (dir_, runnumberstring, "buildings")
    if os.path.isfile(filename) == False:
        print("\n   HDF5 file doesnt exist: %s \n   Skip module\n" % filename)
        return False
    if os.path.exists(output_path_) == False:
        os.mkdir(output_path_)

    path_ = "%s%s" % (output_path_, runnumberstring)

    final_year = 2080

    year_list_hdf5 = []

    print(filename)
    print(os.path.isfile(filename))
    # Access the hdf5 data container
    hdf5_f = h5py.File(filename, "r")
    print(hdf5_f.items())

    bca_names = None
    num_bca = None
    item_names_hdf5 = hdf5_f.items()
    for curr_item_name in item_names_hdf5:
        if curr_item_name[0][:3] == "BC_" and int(curr_item_name[0][-4:]) <= final_year:
            year_list_hdf5.append(int(curr_item_name[0][-4:]))
        elif curr_item_name[0][:18] == "BuildingCategories":
            bca_names = hdf5_f[curr_item_name[0]]["name"]
            num_bca = bca_names.shape[0]

    year_list_hdf5.sort()
    #start_yr = year_list_hdf5[0]
    #final_yr = year_list_hdf5[-1]

    if len(year_list_hdf5) > 1:
        # remove first year
        year_list_hdf5 = year_list_hdf5[1:]

        header = ""
        for yr in year_list_hdf5:
            header += "," + str(yr)
        st0 = time.time()
        for i, yr in enumerate(year_list_hdf5):
            print("Current yr: %i" % yr)
            print("  --> Prepare Data")
            st = time.time()

            key_bc = "BC_%i" % yr
            key_bssh = "BSSH_%i" % yr
            bc = hdf5_f[key_bc][()]
            bssh = hdf5_f[key_bssh][()]
            num_bc = bc.shape[0]  # Number of building classes

            num_bssh = bssh.shape[0]  # Number of building segments
            bssh_bc_idx = bssh["building_classes_index"]

            try:
                # Check if exists
                results_volume_new_buildings * 1
            except BaseException:
                if num_bca is None:
                    num_bca = np.max(bc["building_categories_index"])
                    bca_names = []
                    for n in range(num_bca):
                        bca_names.append("BCA_%i" % (n + 1))
                    bca_names = np.array(bca_names)
                results_volume_new_buildings = np.zeros(
                    (num_bca, len(year_list_hdf5)), dtype="f4"
                )
                results_volume_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_gfa_new_buildings = np.zeros_like(results_volume_new_buildings)
                results_gfa_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_buildings_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_buildings_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_wall_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_wall_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_ceiling_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_ceiling_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_floor_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_floor_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )

                results_volume_internal_structures_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_volume_internal_structures_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )

                results_area_wall_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_wall_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_ceiling_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_ceiling_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_floor_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_floor_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_windows_new_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_windows_renovated_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                results_area_windows_renovated_inc_maint_buildings = np.zeros_like(
                    results_volume_new_buildings
                )

                all_build_results_area_wall_existing_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                all_build_results_area_ceiling_existing_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                all_build_results_area_floor_existing_buildings = np.zeros_like(
                    results_volume_new_buildings
                )
                all_build_results_area_windows_existing_buildings = np.zeros_like(
                    results_volume_new_buildings
                )

            RESULTS_MATRIX = np.zeros((bssh_bc_idx.shape[0], 21))

            # Get Current bc index of buildings
            bc_idx_bc = bc["index"]
            bc_idx_bssh = bc_idx_bc[bssh_bc_idx - 1]

            bca_idx_bssh = bc["building_categories_index"][bc_idx_bssh - 1]

            # Get Former bc index of buildings
            former_bc_idx_bc = bc["prev_index"]
            former_bc_idx_bssh = former_bc_idx_bc[bssh_bc_idx - 1]

            RESULTS_MATRIX[:, 0] = bc_idx_bssh
            RESULTS_MATRIX[:, 1] = former_bc_idx_bssh
            RESULTS_MATRIX[:, 2] = bca_idx_bssh

            # Index of buildings which have been refurbished
            bc_idx_renovation_considered_inc_maint_bc = (
                bc["renovation_facade_year_end"] == yr
            )
            bc_idx_not_maintenance = bc["last_action"] != "maintenance"
            bc_idx_renovation_considered_bc = np.logical_and(
                bc_idx_renovation_considered_inc_maint_bc, bc_idx_not_maintenance
            )

            bc_idx_renovation_considered_inc_maint_bssh = (
                bc_idx_renovation_considered_inc_maint_bc[bssh_bc_idx - 1]
            )
            bc_idx_renovation_considered_bssh = bc_idx_renovation_considered_bc[
                bssh_bc_idx - 1
            ]

            bc_idx_new_buildings_bc = bc["construction_period_end"] == yr
            bc_idx_new_buildings_bssh = bc_idx_new_buildings_bc[bssh_bc_idx - 1]

            ##########################
            #
            # Do the actual Calculation
            #
            ##########################

            U1 = bc["u_value_exterior_walls"][former_bc_idx_bssh - 1]
            U2 = bc["u_value_exterior_walls"][bc_idx_bssh - 1]
            insulation_thickness_wall = np.maximum(
                0,
                np.minimum(
                    max_ins_thickness_wall, ((U1 - U2) / (U1 * U2) * lambda_insulation)
                ),
            )
            # insulation_thickness_wall[insulation_thickness_wall > min_calc_ins_thickness_to_increase] = np.maximum(min_ins_thickness_wall
            #                                                                         , insulation_thickness_wall[insulation_thickness_wall > min_calc_ins_thickness_to_increase])

            idx = np.logical_and(
                insulation_thickness_wall > min_calc_ins_thickness_to_increase,
                insulation_thickness_wall < min_ins_thickness_wall,
            )
            insulation_thickness_wall[idx] = min_ins_thickness_wall
            RESULTS_MATRIX[:, 3] = U1
            RESULTS_MATRIX[:, 4] = np.maximum(min_U_value, U2)
            RESULTS_MATRIX[:, 5] = insulation_thickness_wall

            U1 = bc["u_value_ceiling"][former_bc_idx_bssh - 1]
            U2 = bc["u_value_ceiling"][bc_idx_bssh - 1]
            insulation_thickness_ceiling = np.maximum(
                0,
                np.minimum(
                    max_ins_thickness_ceil, ((U1 - U2) / (U1 * U2) * lambda_insulation)
                ),
            )
            # insulation_thickness_ceiling[insulation_thickness_ceiling > min_calc_ins_thickness_to_increase] = np.maximum(min_ins_thickness_ceil
            #                                                                               , insulation_thickness_ceiling[insulation_thickness_ceiling > min_calc_ins_thickness_to_increase])
            idx = np.logical_and(
                insulation_thickness_ceiling > min_calc_ins_thickness_to_increase,
                insulation_thickness_ceiling < min_ins_thickness_ceil,
            )
            insulation_thickness_ceiling[idx] = min_ins_thickness_ceil
            RESULTS_MATRIX[:, 6] = U1
            RESULTS_MATRIX[:, 7] = np.maximum(min_U_value, U2)
            RESULTS_MATRIX[:, 8] = insulation_thickness_ceiling

            U1 = bc["u_value_floor"][former_bc_idx_bssh - 1]
            U2 = bc["u_value_floor"][bc_idx_bssh - 1]
            insulation_thickness_floor = np.maximum(
                0,
                np.minimum(
                    max_ins_thickness_floor, ((U1 - U2) / (U1 * U2) * lambda_insulation)
                ),
            )
            # insulation_thickness_floor[insulation_thickness_floor > min_calc_ins_thickness_to_increase] = np.maximum(min_ins_thickness_floor
            #                                                                            , insulation_thickness_floor[insulation_thickness_floor > min_calc_ins_thickness_to_increase])
            idx = np.logical_and(
                insulation_thickness_floor > min_calc_ins_thickness_to_increase,
                insulation_thickness_floor < min_ins_thickness_floor,
            )
            insulation_thickness_floor[idx] = min_ins_thickness_floor
            RESULTS_MATRIX[:, 9] = U1
            RESULTS_MATRIX[:, 10] = np.maximum(min_U_value, U2)
            RESULTS_MATRIX[:, 12] = insulation_thickness_floor

            insulation_volume_wall_per_building = non_optimality_factor * (
                vertical_and_horizontal_building_ABGESETZTHEIT
                * bc["aewd"][bc_idx_bssh - 1]
                * insulation_thickness_wall
            )
            insulation_volume_ceiling_per_building = non_optimality_factor * (
                bc["areafloor"][bc_idx_bssh - 1] * insulation_thickness_ceiling
            )
            insulation_volume_floor_per_building = non_optimality_factor * (
                bc["areafloor"][bc_idx_bssh - 1] * insulation_thickness_floor
            )

            insulation_volume_wall_all_buildings_bssh = (
                insulation_volume_wall_per_building * bssh["number_of_buildings"]
            )
            insulation_volume_ceiling_all_buildings_bssh = (
                insulation_volume_ceiling_per_building * bssh["number_of_buildings"]
            )
            insulation_volume_floor_all_buildings_bssh = (
                insulation_volume_floor_per_building * bssh["number_of_buildings"]
            )

            insulation_area_wall_all_buildings_bssh = (
                vertical_and_horizontal_building_ABGESETZTHEIT
                * bc["aewd"][bc_idx_bssh - 1]
                * bssh["number_of_buildings"]
            )
            insulation_area_ceiling_all_buildings_bssh = (
                vertical_and_horizontal_building_ABGESETZTHEIT**0.5
                * bc["areafloor"][bc_idx_bssh - 1]
                * bssh["number_of_buildings"]
            )
            insulation_area_floor_all_buildings_bssh = (
                vertical_and_horizontal_building_ABGESETZTHEIT**0.5
                * bc["areafloor"][bc_idx_bssh - 1]
                * bssh["number_of_buildings"]
            )
            insulation_area_windows_all_buildings_bssh = (
                bc["areawindows"][bc_idx_bssh - 1] * bssh["number_of_buildings"]
            )

            RESULTS_MATRIX[:, 13] = insulation_volume_wall_per_building
            RESULTS_MATRIX[:, 14] = insulation_volume_ceiling_per_building
            RESULTS_MATRIX[:, 15] = insulation_volume_floor_per_building

            gfa_area_all_buildings_bssh = (
                bc["grossfloor_area"][bc_idx_bssh - 1] * bssh["number_of_buildings"]
            )
            number_all_buildings_bssh = bssh["number_of_buildings"]

            # area_buildings_bc = bc['grossfloor_area']

            idx_renovated_all_bca = np.logical_and(
                bc_idx_renovation_considered_bssh == True,
                bc_idx_new_buildings_bssh == False,
            )
            idx_renovated_inc_maint_all_bca = np.logical_and(
                bc_idx_renovation_considered_inc_maint_bssh == True,
                bc_idx_new_buildings_bssh == False,
            )
            idx_new_all_bca = np.logical_and(
                bc_idx_renovation_considered_bssh == True,
                bc_idx_new_buildings_bssh == True,
            )
            idx_existing_all_bca = bc_idx_new_buildings_bssh == False

            RESULTS_MATRIX[idx_renovated_all_bca, 16] = number_all_buildings_bssh[
                idx_renovated_all_bca
            ]
            RESULTS_MATRIX[idx_new_all_bca, 17] = number_all_buildings_bssh[idx_new_all_bca]

            # Estimate Demand of internal structures
            results_volume_internal_walls_per_building = (
                (
                    (bc["length_of_building"] + bc["width_of_building"])
                    * bc["number_of_floors"]
                    * bc["room_height"]
                )
                * thickness_internal_wall
            )[bc_idx_bssh - 1] * bssh["number_of_buildings"]
            results_volume_internal_floors_per_building = (
                bc["grossfloor_area"] * thickness_internal_floor
            )[bc_idx_bssh - 1] * bssh["number_of_buildings"]

            results_volume_internal_structures_renovated_buildings_all_bca = np.zeros(
                (bssh_bc_idx.shape[0])
            )
            results_volume_internal_structures_renovated_buildings_all_bca[
                idx_renovated_all_bca
            ] = (
                results_volume_internal_walls_per_building[idx_renovated_all_bca]
                * share_retrofit_walls
                + results_volume_internal_floors_per_building[idx_renovated_all_bca]
                * share_retrofit_floors
            )
            results_volume_internal_structures_new_buildings_all_bca = np.zeros(
                (bssh_bc_idx.shape[0])
            )
            results_volume_internal_structures_new_buildings_all_bca[
                idx_new_all_bca
            ] = (
                results_volume_internal_walls_per_building[idx_new_all_bca]
                * share_new_buildings_walls
                + results_volume_internal_floors_per_building[idx_new_all_bca]
                * share_new_buildings_floors
            )
            RESULTS_MATRIX[
                :, 19
            ] = results_volume_internal_structures_renovated_buildings_all_bca
            RESULTS_MATRIX[:, 20] = results_volume_internal_structures_new_buildings_all_bca
            ##########################
            #
            # End of actual Calculation
            #
            ##########################

            idx_buildings_some_measure = (RESULTS_MATRIX[:, 16] + RESULTS_MATRIX[:, 17]) > 0
            DUMMY2 = RESULTS_MATRIX[idx_buildings_some_measure, :]
            # print (time.time()-st)
            print("  --> Export Segments")
            output_filename = (
                "%sInsRes_building_segments_%i.csv" % (path_, yr)
            ).replace("\\", "/")
            np.savetxt(output_filename, DUMMY2[:200, :], delimiter=",")  # Export sample dataset
            # print (time.time()-st)
            print("  --> Calc Data per BCA")
            for jj in range(0, num_bca):
                idx_bca = bca_idx_bssh == (jj + 1)
                idx_renovated = np.logical_and(idx_renovated_all_bca, idx_bca)
                idx_renovated_inc_maint = np.logical_and(
                    idx_renovated_inc_maint_all_bca, idx_bca
                )

                idx_new = np.logical_and(idx_new_all_bca, idx_bca)
                idx_exists = np.logical_and(idx_existing_all_bca, idx_bca)

                results_volume_wall_new_buildings[jj, i] = np.sum(
                    insulation_volume_wall_all_buildings_bssh[idx_new]
                )
                results_volume_wall_renovated_buildings[jj, i] = np.sum(
                    insulation_volume_wall_all_buildings_bssh[idx_renovated]
                )
                results_volume_ceiling_new_buildings[jj, i] = np.sum(
                    insulation_volume_ceiling_all_buildings_bssh[idx_new]
                )
                results_volume_ceiling_renovated_buildings[jj, i] = np.sum(
                    insulation_volume_ceiling_all_buildings_bssh[idx_renovated]
                )
                results_volume_floor_new_buildings[jj, i] = np.sum(
                    insulation_volume_floor_all_buildings_bssh[idx_new]
                )
                results_volume_floor_renovated_buildings[jj, i] = np.sum(
                    insulation_volume_floor_all_buildings_bssh[idx_renovated]
                )

                results_gfa_new_buildings[jj, i] = np.sum(
                    gfa_area_all_buildings_bssh[idx_new]
                )
                results_gfa_renovated_buildings[jj, i] = np.sum(
                    gfa_area_all_buildings_bssh[idx_renovated]
                )
                results_buildings_new_buildings[jj, i] = np.sum(
                    number_all_buildings_bssh[idx_new]
                )
                results_buildings_renovated_buildings[jj, i] = np.sum(
                    number_all_buildings_bssh[idx_renovated]
                )

                results_volume_internal_structures_renovated_buildings[jj, i] = np.sum(
                    results_volume_internal_structures_renovated_buildings_all_bca[
                        idx_bca
                    ]
                )
                results_volume_internal_structures_new_buildings[jj, i] = np.sum(
                    results_volume_internal_structures_new_buildings_all_bca[idx_bca]
                )

                results_area_wall_new_buildings[jj, i] = np.sum(
                    insulation_area_wall_all_buildings_bssh[idx_new]
                )
                results_area_wall_renovated_buildings[jj, i] = np.sum(
                    insulation_area_wall_all_buildings_bssh[idx_renovated]
                )
                results_area_ceiling_new_buildings[jj, i] = np.sum(
                    insulation_area_ceiling_all_buildings_bssh[idx_new]
                )
                results_area_ceiling_renovated_buildings[jj, i] = np.sum(
                    insulation_area_ceiling_all_buildings_bssh[idx_renovated]
                )
                results_area_floor_new_buildings[jj, i] = np.sum(
                    insulation_area_floor_all_buildings_bssh[idx_new]
                )
                results_area_floor_renovated_buildings[jj, i] = np.sum(
                    insulation_area_floor_all_buildings_bssh[idx_renovated]
                )
                results_area_windows_new_buildings[jj, i] = np.sum(
                    insulation_area_windows_all_buildings_bssh[idx_new]
                )
                results_area_windows_renovated_buildings[jj, i] = np.sum(
                    insulation_area_windows_all_buildings_bssh[idx_renovated]
                )

                results_area_windows_renovated_inc_maint_buildings[jj, i] = np.sum(
                    insulation_area_windows_all_buildings_bssh[idx_renovated_inc_maint]
                )

                all_build_results_area_wall_existing_buildings[jj, i] = np.sum(
                    insulation_area_wall_all_buildings_bssh[idx_exists]
                )
                all_build_results_area_ceiling_existing_buildings[jj, i] = np.sum(
                    insulation_area_ceiling_all_buildings_bssh[idx_exists]
                )
                all_build_results_area_floor_existing_buildings[jj, i] = np.sum(
                    insulation_area_floor_all_buildings_bssh[idx_exists]
                )
                all_build_results_area_windows_existing_buildings[jj, i] = np.sum(
                    insulation_area_windows_all_buildings_bssh[idx_exists]
                )

            print(time.time() - st)
            print("NExt Year")
        print("Total_cal_time:")
        print(time.time() - st0)

        results_volume_new_buildings[:, :] = (
            results_volume_wall_new_buildings[:, :]
            + results_volume_ceiling_new_buildings[:, :]
            + results_volume_floor_new_buildings[:, :]
        )

        results_volume_renovated_buildings[:, :] = (
            results_volume_ceiling_renovated_buildings[:, :]
            + results_volume_ceiling_renovated_buildings[:, :]
            + results_volume_floor_renovated_buildings[:, :]
        )

        print("EXPORTING CSV")
        if 1 == 1:
            path_ = path_.replace("\\", "/")
            output_filename = (
                "%sInsRes_external_insulation_volume_wall_new_buildings.csv" % (path_)
            )
            print(output_filename)
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_wall_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_wall_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_wall_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_ceiling_new_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_ceiling_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_ceiling_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, results_volume_ceiling_renovated_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_floor_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_floor_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_floor_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_floor_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_renovated_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_volume_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = "%sInsRes_gfa_new_buildings.csv" % (path_)
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_gfa_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = "%sInsRes_gfa_renovated_buildings.csv" % (path_)
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_gfa_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = "%sInsRes_num_buildings_new_buildings.csv" % (path_)
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_buildings_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = "%sInsRes_num_buildings_renovated_buildings.csv" % (path_)
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_buildings_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_all_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_new_buildings
                        + results_volume_renovated_buildings,
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_all_buildings_single_line.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        ["All BCA"],
                        np.expand_dims(
                            np.sum(
                                results_volume_new_buildings
                                + results_volume_renovated_buildings,
                                axis=0,
                            ),
                            axis=0,
                        ),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_per_gfa_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_renovated_buildings
                        / np.maximum(0.001, results_gfa_renovated_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_per_gfa_new_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_new_buildings
                        / np.maximum(0.001, results_gfa_new_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_volume_per_gfa_all_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        (
                            results_volume_new_buildings
                            + results_volume_renovated_buildings
                        )
                        / np.maximum(
                            0.001,
                            results_gfa_new_buildings + results_gfa_renovated_buildings,
                        ),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, results_volume_internal_structures_new_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_renovated_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, results_volume_internal_structures_renovated_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_all_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_internal_structures_new_buildings
                        + results_volume_internal_structures_renovated_buildings,
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_per_gfa_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_internal_structures_renovated_buildings
                        / np.maximum(0.001, results_gfa_renovated_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_per_gfa_new_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_internal_structures_new_buildings
                        / np.maximum(0.001, results_gfa_new_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_internal_insulation_volume_per_gfa_all_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        (
                            results_volume_internal_structures_new_buildings
                            + results_volume_internal_structures_renovated_buildings
                        )
                        / np.maximum(
                            0.001,
                            results_gfa_new_buildings + results_gfa_renovated_buildings,
                        ),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            # results_volume_new_buildings + results_volume_renovated_buildings

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_new_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_new_buildings
                        + results_volume_internal_structures_new_buildings,
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_renovated_buildings
                        + results_volume_internal_structures_renovated_buildings,
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_all_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        results_volume_new_buildings
                        + results_volume_internal_structures_new_buildings
                        + results_volume_renovated_buildings
                        + results_volume_internal_structures_renovated_buildings,
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_per_gfa_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        (
                            results_volume_renovated_buildings
                            + results_volume_internal_structures_renovated_buildings
                        )
                        / np.maximum(0.001, results_gfa_renovated_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_per_gfa_new_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        (
                            results_volume_new_buildings
                            + results_volume_internal_structures_new_buildings
                        )
                        / np.maximum(0.001, results_gfa_new_buildings),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_per_gfa_all_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        bca_names,
                        (
                            results_volume_new_buildings
                            + results_volume_internal_structures_new_buildings
                            + results_volume_renovated_buildings
                            + results_volume_internal_structures_renovated_buildings
                        )
                        / np.maximum(
                            0.001,
                            results_gfa_new_buildings + results_gfa_renovated_buildings,
                        ),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_AND_internal_insulation_volume_all_buildings_single_line.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (
                        ["All BCA"],
                        np.expand_dims(
                            np.sum(
                                results_volume_new_buildings
                                + results_volume_renovated_buildings,
                                axis=0,
                            )
                            + np.sum(
                                results_volume_internal_structures_new_buildings
                                + results_volume_internal_structures_renovated_buildings,
                                axis=0,
                            ),
                            axis=0,
                        ),
                    )
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_area_wall_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_wall_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_wall_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_wall_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_wall_all_existing_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, all_build_results_area_wall_existing_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_area_ceiling_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_ceiling_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_ceiling_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_ceiling_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_ceiling_all_existing_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, all_build_results_area_ceiling_existing_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = (
                "%sInsRes_external_insulation_area_floor_new_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_floor_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_floor_renovated_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_floor_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_insulation_area_floor_all_existing_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, all_build_results_area_floor_existing_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

            output_filename = "%sInsRes_external_windows_area_new_buildings.csv" % (
                path_
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_windows_new_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_windows_area_renovated_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack((bca_names, results_area_windows_renovated_buildings)),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_windows_area_all_existing_buildings.csv" % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, all_build_results_area_windows_existing_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )
            output_filename = (
                "%sInsRes_external_windows_area_renovated_incMaint_buildings.csv"
                % (path_)
            )
            np.savetxt(
                output_filename,
                np.column_stack(
                    (bca_names, results_area_windows_renovated_inc_maint_buildings)
                ),
                fmt="%s",
                delimiter=",",
                header=header,
            )

    hdf5_f.close()
    print("  -----  THE END  --------")
    print("\n\n#####################")
    print("#####################\n\n")


if __name__ == "__main__":

    path_ = "./testdata/001_"
    path_ = path_.replace("\\", "/")

    SCENARIO_INFO = {}
    CreateInsulationData(SCENARIO_INFO, path_)
