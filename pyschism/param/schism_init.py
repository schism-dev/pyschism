import pathlib
import re
import tempfile
import urllib.request

import f90nml


def typecast_value(value: str):
    if "real(" in value:
        value, exponent = value.split("(")[-1].split(",")[0].lower().split("d")
        return float(f"{value}e{exponent}")
    elif "_rkind" in value:
        return float(value.split("_rkind")[0])
    elif "'" in value or '"' in value:
        return value[1:-1]
    elif re.match("^[0-9].[0-9]*[de][+-]*[0-9]*", value):
        value, exponent = value.split("d")
        return float(f"{value}e{exponent}")
    elif "pi" in value:
        return value
    elif re.match("[+-]?([0-9]*[.])?[0-9]+", value):
        if "." in value:
            return float(value)
        else:
            return int(value)
    else:
        return value


class NamelistParser:
    def __new__(cls, schism_init_f90: str, component):
        schism_init_f90 = schism_init_f90.split("\n")
        for i, line in enumerate(schism_init_f90):
            line = line.lower()
            if f"namelist/{component}/" in line.replace(" ", ""):
                declared_variables = (
                    line.strip()
                    .split("namelist")[1]
                    .strip()
                    .split(f"/{component}/")[1]
                    .strip()
                    .split("&")[0]
                    .strip()
                    .split(",")[:-1]
                )
                break
        line = "&&"

        while line[0] == line[-1] == "&":
            line = schism_init_f90[i + 1].strip()
            declared_variables.extend(
                [var for var in line.replace("&", "").strip().split(",") if var]
            )
            i += 1
        defaults = {var: None for var in declared_variables}
        for var in declared_variables:
            for line in schism_init_f90[i + 1 :]:
                if "!" in line:
                    continue
                for source_var in line.split(";"):
                    source_var = source_var.strip()
                    if (
                        re.match(
                            rf"^{var}[ ]*=[ ]*[-_A-Za-z0-9 .\(\),]*[;]*", source_var
                        )
                        is not None
                    ):
                        value = (
                            [x for x in line.replace(" ", "").split(";") if var in x]
                            .pop()
                            .split("=")[-1]
                        )
                        defaults.update({var: typecast_value(value)})
                        break

        return defaults


class CoreParser(NamelistParser):
    def __new__(cls, schism_init_f90):
        return super().__new__(cls, schism_init_f90, "core")


class OptParser(NamelistParser):
    def __new__(cls, schism_init_f90):
        return super().__new__(cls, schism_init_f90, "opt")


class SchoutParser(NamelistParser):
    def __new__(cls, schism_init_f90):
        return super().__new__(cls, schism_init_f90, "schout")


class SchismInitParser:
    def __init__(self, branch="master"):
        url = f"https://raw.githubusercontent.com/schism-dev/schism/{branch}/src/Hydro/schism_init.F90"
        response = urllib.request.urlopen(url)
        self.schism_init_f90 = response.read().decode("utf-8")

    @property
    def opt(self):
        return OptParser(self.schism_init_f90)

    @property
    def core(self):
        return CoreParser(self.schism_init_f90)

    @property
    def schout(self):
        return SchoutParser(self.schism_init_f90)


def patch_sample_param_to_avoid_f90nml_bug(schism_param_sample):
    pattern = r"^[! _A-Za-z0-9]+\([*\d]*(\d+)[^\d]*\)[ ]*="
    var_collection = []
    for i, line in enumerate(schism_param_sample):
        if re.match(pattern, line):
            if line[0] != "!":
                schism_param_sample[i] = line = "!" + line
            varname = line.strip().split("=")[0].strip()[1:].strip().split("(")[0]
            var_collection.append(varname)
    declared_vars = list(set(var_collection))
    for var in declared_vars:
        varcount = var_collection.count(var)
        for i, line in enumerate(list(schism_param_sample)):
            if f"{var}({varcount})" in line:
                if varcount > 1:
                    schism_param_sample.insert(
                        i + 1, f"{var}(1:{varcount}) = " + ", ".join(varcount * ["0"])
                    )
                else:
                    schism_param_sample.insert(i + 1, f"{var}({varcount}) = 0")
                break


class Singleton(type):
    def __init__(cls, name, bases, dict):
        super(Singleton, cls).__init__(name, bases, dict)
        cls.instance = None

    def __call__(cls, *args, **kw):
        if cls.instance is None:
            cls.instance = super(Singleton, cls).__call__(*args, **kw)
        return cls.instance


class GitParamTemplate(metaclass=Singleton):
    def __init__(self, branch="master"):
        url = f"https://raw.githubusercontent.com/schism-dev/schism/{branch}/sample_inputs/param.nml"
        response = urllib.request.urlopen(url)
        schism_param_sample = response.read().decode("utf-8").split("\n")
        patch_sample_param_to_avoid_f90nml_bug(schism_param_sample)
        tmpfile1 = tempfile.NamedTemporaryFile()
        with open(tmpfile1.name, "w") as f:
            f.write("\n".join(schism_param_sample))
        base_param = f90nml.read(tmpfile1.name)

        def update_param_component(name, component):
            for source_code_var, source_code_value in component.items():
                if source_code_var in base_param[name]:
                    if isinstance(base_param[name][source_code_var], list):
                        continue
                    if source_code_value != base_param[name][source_code_var]:
                        base_param[name][source_code_var] = source_code_value

        schism_init_f90 = SchismInitParser()
        update_param_component("core", schism_init_f90.core)
        update_param_component("opt", schism_init_f90.opt)
        update_param_component("schout", schism_init_f90.schout)
        tmpfile2 = tempfile.NamedTemporaryFile()
        f90nml.patch(tmpfile1.name, base_param, tmpfile2.name)
        with open(tmpfile2.name) as f:
            self.schism_param_sample = "".join(f.readlines())
        self.tmpfile = tmpfile2
        self.schism_init_f90 = schism_init_f90

    def __str__(self):
        return self.schism_param_sample

    def write(self, path):
        pathlib.Path(path).touch()
        f90nml.patch(self.path, str(self), path)

    @property
    def core(self):
        return self.schism_init_f90.core

    @property
    def opt(self):
        return self.schism_init_f90.opt

    @property
    def schout(self):
        return self.schism_init_f90.schout

    @property
    def path(self):
        return pathlib.Path(self.tmpfile.name)
