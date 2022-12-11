import argparse
import datetime
import os
import re
import yaml
import pprint


def load_argumets():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tpr_file',
                        type=str,
                        help='Tpr file from gromacs with metadata.')
    parser.add_argument('--gro_file',
                        type=str,
                        help='Gro file from gromacs with metadata.')
    parser.add_argument('--cpt_file',
                        type=str,
                        help='cpt file from gromacs with metadata.',
                        default=None)
    parser.add_argument('--print_metadata',
                        type=str,
                        choices=("json", "yaml"),
                        default="yaml",
                        help="Print extracted metadata as json or yaml.")
    args = parser.parse_args()

    if not os.path.isfile(args.tpr_file):
        exit(f"\nERROR! There is no file {args.tpr_file}!\n")
    if not os.path.isfile(args.gro_file):
        exit(f"\nERROR! There is no file {args.gro_file}!\n")

    return args


def load_tpr_lines(tpr_file: str) -> list:
    os.system(f"gmx dump -s {tpr_file} 1> metadata_gmx.tmp")
    tpr_lines = [line.split() for line in open("metadata_gmx.tmp", "r").read().lower().splitlines()]
    #os.system("rm metadata_gmx.tmp")
    return tpr_lines


def extract_main_information(args) -> dict:
    def _list_of_molecules():
        shortcuts = {'ala': 'A',
                     'arg': 'R',
                     'asn': 'N',
                     'asp': 'D',
                     'asx': 'B',
                     'cys': 'C',
                     'glu': 'E',
                     'gln': 'Q',
                     'glx': 'Z',
                     'gly': 'G',
                     'his': 'H',
                     'ile': 'I',
                     'leu': 'L',
                     'lys': 'K',
                     'met': 'M',
                     'phe': 'F',
                     'pro': 'P',
                     'ser': 'S',
                     'thr': 'T',
                     'trp': 'W',
                     'tyr': 'Y',
                     'val': 'V'}
        mols_data = {}
        mcounter = 0
        for i, sl in enumerate(tpr_lines):
            if sl == ["molblock", f"({mcounter}):"]:
                mcounter += 1
                mname = tpr_lines[i + 1][-1][1:-1]
                nmols = tpr_lines[i + 2][-1]
                mols_data[f"molecule {mcounter}"] = {"name": mname,
                                                     "count": nmols}
        for x in range(mcounter):
            for index, sl in enumerate(tpr_lines):
                if bool(re.search("moltype \((\d+)\):", " ".join(sl))) and tpr_lines[index + 1] == [f"name=\"{mols_data[f'molecule {x + 1}']['name']}\""]:
                    while not bool(re.search("residue \((\d+)\)", " ".join(tpr_lines[index]))):
                        index += 1
                    index += 1
                    residues = []
                    while tpr_lines[index][0][:7] == "residue":
                        residues.append("".join(tpr_lines[index]).split("=")[2][1:-4])
                        index += 1
                    if len(list(set(shortcuts.keys()).intersection(set(residues)))):
                        mols_data[f"molecule {x + 1}"]["residues"] = [shortcuts[res] for res in residues]
        main_information["molecules"] = mols_data

    main_information = {"force field": "probably has to be set by the user"}

    os.system(f"gmx dump -s {args.tpr_file} 2>&1 | grep \"VERSION\" > version_gmx.tmp")
    main_information["gromacs version"] = open("version_gmx.tmp", "r").read().splitlines()[0].split()[4]
    os.system("rm version_gmx.tmp")
    os.system(f"gmx dump -s {args.tpr_file} 2>&1 > version_gmx.tmp")
    main_information["size and shape of the simulation box"] = [float(n) for n in
                                                                open(args.gro_file, "r").readlines()[
                                                                    -1].rstrip().split()]
    os.system("rm version_gmx.tmp")
    integrator = get_value_for_key("integrator")
    if integrator in ["steep", "cg"]:
        main_information["type of simulation"] = "energy minimization"
    elif integrator == "sd" or integrator[:2] == "md":
        main_information["type of simulation"] = "molecular dynamics"
        main_information["simulation time step [ps]"] = get_value_for_key("dt")
        if args.cpt_file is None:
            main_information["simulation length [ps]"] = main_information[
                                                             "simulation time step [ps]"] * get_value_for_key("nsteps")
        else:
            pass
            os.system(f"gmx check -f {args.cpt_file} 2>&1 | grep \"Last frame\" > frame.tmp")
            main_information["simulation length [ps]"] = float(open("frame.tmp", "r").read().splitlines()[1].split()[4])
            os.system("rm frame.tmp")
    tcoupl = get_value_for_key("tcoupl")
    pcoupl = get_value_for_key("pcoupl")
    if tcoupl == pcoupl == "no":
        main_information["statistical ensamble"] = "NVE (microcanonical)"
    elif tcoupl in ["nose-hoover", "v-rescale", "berendsen", "andersen", "andersen-massive"] and pcoupl == "no":
        main_information["statistical ensamble"] = "NVT (canonical)"
        main_information["reference temperature"] = find_all_values_for_key("ref-t")
    elif tcoupl in ["nose-hoover", "v-rescale", "berendsen", "andersen", "andersen-massive"] and pcoupl in ["berendsen", "parrinello-rahman"]:
        main_information["statistical ensamble"] = "NpT (isothermal-isobaric)"
        main_information["reference temperature"] = find_all_values_for_key("ref-t")
        main_information["reference pressure"] = find_matrix_with_key("ref-p")
    _list_of_molecules()
    store_values_for_keys(target=main_information,
                          keys=(("free energy calculation", "free-energy"),
                                ("umbrella sampling", "pull"),
                                ("AWH adaptive biasing", "awh")))
    return main_information, tcoupl, pcoupl


def extract_published_information(tcoupl: str,
                                  pcoupl: str) -> dict:
    def _find_tc_grps():
        for sl in tpr_lines:
            if sl[0] == "grp[t-coupling":
                nr = sl[2][3:-1]
                name = " ".join(sl[4:])[:-1]
                published_information["thermostat"]["tc-grps"] = {"nr": nr, "name": name}
                break
        else:
            print("Warning, no value for tc-grps founded!")

    published_information = {}
    if main_information["type of simulation"] == "energy minimization":
        store_values_for_keys(target=published_information,
                              keys=("emtol",
                                    "emstep"))
    if tcoupl != "no":
        store_values_for_keys(target=published_information,
                              keys=("tcoupl",
                                    "nsttcouple"),
                              group="thermostat")
        published_information["thermostat"]["tau-t"] = find_all_values_for_key("tau-t")
        _find_tc_grps()
    if pcoupl != "no":
        store_values_for_keys(target=published_information,
                              keys=("pcoupl",
                                    "refcoord-scaling",
                                    "pcoupltype",
                                    "tau-p"),
                              group="barostat")
        published_information["barostat"]["compressibility"] = find_matrix_with_key("compressibility")
    store_values_for_keys(target=published_information,
                          keys=("rvdw",
                                "vdw-type",
                                "rvdw-switch",
                                "vdw-modifier",
                                "dispcorr"),
                          group="van der Waals interactions")
    store_values_for_keys(target=published_information,
                          keys=("coulombtype",
                                "coulomb-modifier",
                                "rcoulomb",
                                "epsilon-r",
                                "epsilon-rf"),
                          group="electrostatic interactions")
    store_values_for_keys(target=published_information,
                          keys=("cutoff-scheme",
                                "nstlist",
                                "pbc",
                                "rlist"),
                          group="neighbour list")
    return published_information


def extract_detailed_information() -> dict:
    detailed_information = {}
    store_values_for_keys(target=detailed_information,
                          keys=("constraint-algorithm",
                                "constraints",
                                "lincs-order",
                                "lincs-iter",
                                "comm-grps",
                                "nstcomm",
                                "comm-mode"))
    if published_information["electrostatic interactions"]["coulombtype"] == "pme":
        store_values_for_keys(target=detailed_information,
                              keys=("fourierspacing",))
    if main_information["free energy calculation"] == "yes":
        store_values_for_keys(target=detailed_information,
                              keys=("init-lambda",
                                    "delta-lambda",
                                    "sc-alpha",
                                    "sc-power",
                                    "sc-sigma"))
    return detailed_information


def add_umbrella_sampling_information():
    def _find_dim(index: int,
                  i: int):
        for index, sl in enumerate(tpr_lines[index:], start=index):
            if bool(re.search("dim \((\d+)\):", " ".join(sl))):
                dim = []
                for n in range(1, int(sl[1][1:-2]) + 1):
                    dim.append(convert_str_to_int_float("".join(tpr_lines[index + n]).split("=")[1]))
                break
        published_information["umbrella sampling"][f"pull-coord{i}-dim"] = dim

    def _find_vec(index: int,
                  i: int):
        for index, sl in enumerate(tpr_lines[index:], start=index):
            if bool(re.search("vec \((\d+)\):", " ".join(sl))):
                vec = []
                for n in range(1, int(sl[1][1:-2]) + 1):
                    vec.append(convert_str_to_int_float("".join(tpr_lines[index + n]).split("=")[1]))
                break
        published_information["umbrella sampling"][f"pull-coord{i}-vec"] = vec

    def _find_groups(index: int,
                     i: int):
        while not bool(re.search("group\[(\d+)\]", tpr_lines[index][0])):
            index += 1
        groups = []
        while bool(re.search("group\[(\d+)\]", tpr_lines[index][0])):
            groups.append(convert_str_to_int_float(tpr_lines[index][2]))
            index += 1
        published_information["umbrella sampling"][f"pull-coord{i}-groups"] = groups

    published_information["umbrella sampling"] = {}
    store_values_for_keys(target=published_information["umbrella sampling"],
                          keys=("pull-ncoords",
                                "pull-ngroups"))
    for i in range(published_information["umbrella sampling"]["pull-ncoords"]):
        for index, sl in enumerate(tpr_lines):
            if sl == ["pull-coord", f"{i}:"]:
                _find_dim(index, i)
                _find_vec(index, i)
                _find_groups(index, i)

                store_first_values_for_keys_after_index(target=detailed_information,
                                                        names_keys=((f"pull-coord{i}-start", "start"),
                                                                    (f"pull-coord{i}-init", "init")),
                                                        index=index)
                store_first_values_for_keys_after_index(target=published_information["umbrella sampling"],
                                                        names_keys=((f"pull-coord{i}-type", "type"),
                                                                    (f"pull-coord{i}-geometry", "geometry"),
                                                                    (f"pull-coord{i}-rate", "rate"),
                                                                    (f"pull-coord{i}-k", "k")),
                                                        index=index)
                if published_information["umbrella sampling"][f"pull-coord0-geometry"] == "cylinder":
                    store_values_for_keys(target=published_information["umbrella sampling"],
                                          keys=("pull-cylinder-r",))


def add_awh_information():
    published_information["AWH adaptive biasing"] = {}
    store_values_for_keys(target=published_information["AWH adaptive biasing"],
                          keys=("awh-nbias",))
    store_values_for_keys(target=detailed_information,
                          keys=("awh-potential",
                                "awh-share-bias-multisim"))
    for awh_i in range(1, published_information["AWH adaptive biasing"]["awh-nbias"] + 1):
        store_values_for_keys(target=published_information["AWH adaptive biasing"],
                              keys=(f"awh{awh_i}-ndim",
                                    f"awh{awh_i}-error-init"))
        store_values_for_keys(target=detailed_information,
                              keys=(f"awh{awh_i}-target",
                                    f"awh{awh_i}-growth",
                                    f"awh{awh_i}-equilibrate-histogram"))

        for dim_i in range(1, published_information["AWH adaptive biasing"][f"awh{awh_i}-ndim"] + 1):
            for index, sl in enumerate(tpr_lines):
                if sl == [f"awh{awh_i}-dim{dim_i}:"]:
                    store_first_values_for_keys_after_index(target=published_information["AWH adaptive biasing"],
                                                            names_keys=((f"awh{awh_i}-dim{dim_i}-cover-diameter", "cover-diameter"),
                                                                        (f"awh{awh_i}-dim{dim_i}-diffusion", "diffusion"),
                                                                        (f"awh{awh_i}-dim{dim_i}-coord-index", "coord-index"),
                                                                        (f"awh{awh_i}-dim{dim_i}-start", "start"),
                                                                        (f"awh{awh_i}-dim{dim_i}-end", "end"),
                                                                        (f"awh{awh_i}-dim{dim_i}-force-constant", "force-constant")),
                                                            index=index)


def convert_str_to_int_float(value):
    try:
        value = int(value)
    except ValueError:
        try:
            value = float(value)
        except ValueError:
            pass
    return value


def get_value_for_key(key: str):
    filtered = [sl for sl in tpr_lines if [key, "="] == sl[:2] and len(sl) == 3]
    if len(filtered) > 1:
        print(f"Warning! More lines with keyword {key}!")
        return False
    elif len(filtered) == 0:
        print(f"Warning! No value for keyword {key}!")
        return False
    else:
        value = filtered[0][2]
        converted_value = convert_str_to_int_float(value)
        return converted_value


def store_values_for_keys(target: dict,
                          keys: tuple,
                          group: str = None):
    """If the key is a string, the value is stored under the same name as in the tpr file.
    If the key is tuple, the first element is used as the name to store in the metadata
    and the value in the tpr file is looked up by the second element."""
    if group:
        target[group] = {}
        target = target[group]
    for key in keys:
        if isinstance(key, str):
            name = key
        elif isinstance(key, tuple):
            name = key[0]
            key = key[1]
        filtered = [sl for sl in tpr_lines if [key, "="] == sl[:2] and len(sl) == 3]
        if len(filtered) > 1:
            print(f"Warning! More lines with keyword {key}!")
        elif len(filtered) == 0:
            print(f"Warning! No value for keyword {key}!")
        else:
            value = filtered[0][2]
            converted_value = convert_str_to_int_float(value)
            target[name] = converted_value


def store_first_values_for_keys_after_index(target: dict,
                                            names_keys: tuple,
                                            index: int):
    for name, key in names_keys:
        for sl in tpr_lines[index:]:
            if sl[0] == key:
                target[name] = convert_str_to_int_float(sl[2])
                break
        else:
            print(f"Warning! No value for keyword {key}!")


def find_all_values_for_key(key: str) -> list:
    for sl in tpr_lines:
        if sl[0] == f"{key}:":
            return [convert_str_to_int_float(value) for value in sl[1:]]
    print(f"Warning! No values for keyword {key}")


def find_matrix_with_key(key: str) -> list:
    filtered = [sl for sl in tpr_lines if sl[0][:len(key)] == key]
    if len(filtered) != 4:
        exit(f"Wrong number of lines for keyword {key}")
    matrix = [[], [], []]
    for line_count in range(1, 4):
        for index in range(2, 5):
            matrix[line_count - 1].append(float(filtered[line_count][index].replace("}", "").replace(",", "")))
    return matrix


def add_additional_info():
    metadata["_version"] = "1.0.0"
    metadata["_created"] = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")


def print_metadata(json, style):
    if style == "json":
        pprint(json)
    elif style == "yaml":
        print(yaml.dump(json,
                        default_flow_style=False,
                        indent=4,
                        allow_unicode=True)
              .replace("'", ""))


if __name__ == "__main__":
    args = load_argumets()
    tpr_lines = load_tpr_lines(args.tpr_file)

    #exit()
    main_information, tcoupl, pcoupl = extract_main_information(args)
    published_information = extract_published_information(tcoupl, pcoupl)
    detailed_information = extract_detailed_information()

    if main_information["umbrella sampling"] not in ["no", "false"]:
        add_umbrella_sampling_information()

    if main_information["AWH adaptive biasing"] not in ["no", "false"]:
        add_awh_information()

    metadata = {"main_information": main_information,
                "published_information": published_information,
                "detailed_information": detailed_information}

    add_additional_info()
    print_metadata(metadata, args.print_metadata)
