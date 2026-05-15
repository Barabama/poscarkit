"""Form builders — each returns (frame, get_args_callable)."""

import argparse
import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path


# ------------------------------------------------------------------ #
#  Shared widget helpers                                             #
# ------------------------------------------------------------------ #

def _cfg_get(cfg: dict, key: str, default=None):
    return cfg.get(key, default)


def _row(parent, label, widget, **pack_kw):
    """Place a label+widget row in a frame."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2, **pack_kw)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    widget.pack(side=tk.LEFT, fill=tk.X, expand=True)
    return f


def _entry(parent, cfg, key, default=""):
    """Create a labelled Entry with value from config or default."""
    var = tk.StringVar(value=str(_cfg_get(cfg, key, default)))
    w = tk.Entry(parent, textvariable=var)
    return w, var


def _file_row(parent, label, var, cfg, key):
    """Entry + Browse button row."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    e = tk.Entry(f, textvariable=var)
    e.pack(side=tk.LEFT, fill=tk.X, expand=True)
    tk.Button(
        f, text="Browse", command=lambda: _browse_file(var)
    ).pack(side=tk.LEFT, padx=(4, 0))
    return f


def _dir_row(parent, label, var):
    """Entry + Browse folder button row."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    e = tk.Entry(f, textvariable=var)
    e.pack(side=tk.LEFT, fill=tk.X, expand=True)
    tk.Button(
        f, text="Browse", command=lambda: _browse_dir(var)
    ).pack(side=tk.LEFT, padx=(4, 0))
    return f


def _browse_file(var: tk.StringVar):
    path = filedialog.askopenfilename(
        title="Select file",
        filetypes=[("VASP/CSV/XLSX files", "*.vasp *.csv *.xlsx"),
                   ("All files", "*.*")],
    )
    if path:
        var.set(path)


def _browse_dir(var: tk.StringVar):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        var.set(path)


def _checkbox(parent, label, var, default=False):
    """Create a checkbox."""
    v = tk.BooleanVar(value=default)
    tk.Checkbutton(parent, text=label, variable=v).pack(anchor="w")
    return v


def _combo(parent, label, values, cfg, key, default=""):
    """Labelled Combobox."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    var = tk.StringVar(value=str(_cfg_get(cfg, key, default)))
    cb = ttk.Combobox(f, textvariable=var, values=values, state="readonly")
    cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
    return cb, var


def _int_entries(parent, label, count, cfg, key, defaults):
    """Row of N integer entries."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    vals = _cfg_get(cfg, key, list(defaults))
    vars_ = []
    for i in range(count):
        v = tk.StringVar(value=str(vals[i] if i < len(vals) else defaults[i]))
        tk.Entry(f, textvariable=v, width=5).pack(side=tk.LEFT, padx=1)
        vars_.append(v)
    return vars_


# ------------------------------------------------------------------ #
#  Modeling                                                          #
# ------------------------------------------------------------------ #

def modeling_form(parent, cfg: dict):
    name_w, name_var = _entry(parent, cfg, "name", "modeling")
    poscar_var = tk.StringVar(value=str(_cfg_get(cfg, "poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var, cfg, "poscar")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    _, phase_var = _combo(parent, "Phase", ["", "FCC", "BCC", "HCP"], cfg, "phase")
    config_var = tk.StringVar(value="config.toml")
    _file_row(parent, "Config file", config_var, cfg, "config")
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    _, seeds_var = _entry(parent, cfg, "shuffle_seeds", "42")
    batch_w, batch_var = _entry(parent, cfg, "batch_size", "1")
    sqs_var = _checkbox(parent, "Enable SQS", _cfg_get(cfg, "enable_sqs", False))
    iter_w, iter_var = _entry(parent, cfg, "iterations", "10000000")

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            phase=phase_var.get(),
            config=config_var.get() or None,
            factors=[int(v.get()) for v in factors_vars],
            seeds=[int(s.strip()) for s in seeds_var.get().split(",") if s.strip()] or [None],
            batch_size=int(batch_var.get() or "1"),
            enable_sqs=sqs_var.get(),
            iterations=int(float(iter_var.get() or "1e7")),
        )

    # Build the actual form frame
    return get_args


# ------------------------------------------------------------------ #
#  Import SOFs                                                       #
# ------------------------------------------------------------------ #

def import_sofs_form(parent, cfg: dict):
    csv_var = tk.StringVar()
    _file_row(parent, "CSV/XLSX file", csv_var, cfg, "")
    _, phase_var = _combo(parent, "Phase", ["", "FCC", "BCC", "HCP"], cfg, "phase")
    temps_var = tk.StringVar()
    f_t = tk.Frame(parent)
    f_t.pack(fill=tk.X, pady=2)
    tk.Label(f_t, text="Temperatures", width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    tk.Entry(f_t, textvariable=temps_var).pack(side=tk.LEFT, fill=tk.X, expand=True)
    tk.Label(f_t, text="(space-separated)", fg="gray").pack(side=tk.LEFT, padx=(4, 0))

    map_var = tk.StringVar()
    f_m = tk.Frame(parent)
    f_m.pack(fill=tk.X, pady=2)
    tk.Label(f_m, text="Sublattice map", width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    tk.Entry(f_m, textvariable=map_var).pack(side=tk.LEFT, fill=tk.X, expand=True)
    tk.Label(f_m, text="e.g. 1:1a,2:3c", fg="gray").pack(side=tk.LEFT, padx=(4, 0))

    outdir_var = tk.StringVar(value="output")
    _dir_row(parent, "Output dir", outdir_var)
    name_w, name_var = _entry(parent, cfg, "name", "modeling")
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))

    def get_args():
        temps = [float(t) for t in temps_var.get().split()] if temps_var.get().strip() else []
        return argparse.Namespace(
            csv=csv_var.get(),
            config="config.toml",
            phase=phase_var.get(),
            temperatures=temps,
            sublattice_map=map_var.get() or None,
            outdir=outdir_var.get(),
            name=name_var.get(),
            factors=[int(v.get()) for v in factors_vars],
            output="run",
            seeds=None,
            batch_size=1,
            enable_sqs=False,
            iterations=int(1e7),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Count CN                                                          #
# ------------------------------------------------------------------ #

def countcn_form(parent, cfg: dict):
    name_w, name_var = _entry(parent, cfg, "name", "countcn")
    poscar_var = tk.StringVar(value=str(_cfg_get(cfg, "poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var, cfg, "poscar")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    mult_w, mult_var = _entry(parent, cfg, "cutoff_mult", "1.1")
    par_w, par_var = _entry(parent, cfg, "parallel", "2")
    ase_var = _checkbox(parent, "Use ASE backend", _cfg_get(cfg, "by_ase", False))
    pbc_var = _checkbox(parent, "Periodic boundary (PBC)", _cfg_get(cfg, "pbc", False))

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            cutoff_mult=float(mult_var.get() or "1.1"),
            parallel=int(par_var.get() or "2"),
            by_ase=ase_var.get(),
            pbc=pbc_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Slice                                                             #
# ------------------------------------------------------------------ #

def slice_form(parent, cfg: dict):
    name_w, name_var = _entry(parent, cfg, "name", "slice")
    poscar_var = tk.StringVar(value=str(_cfg_get(cfg, "poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var, cfg, "poscar")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    miller_vars = _int_entries(parent, "Miller index", 3, cfg,
                               "slice_direction", (0, 0, 1))

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            miller_index=[int(v.get()) for v in miller_vars],
        )

    return get_args


# ------------------------------------------------------------------ #
#  Slice + CN                                                        #
# ------------------------------------------------------------------ #

def slice_to_countcn_form(parent, cfg: dict):
    name_w, name_var = _entry(parent, cfg, "name", "slice-to-countcn")
    poscar_var = tk.StringVar(value=str(_cfg_get(cfg, "poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var, cfg, "poscar")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    miller_vars = _int_entries(parent, "Miller index", 3, cfg,
                               "slice_direction", (0, 0, 1))
    pbc_var = _checkbox(parent, "Periodic boundary (PBC)", _cfg_get(cfg, "pbc", False))

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            miller_index=[int(v.get()) for v in miller_vars],
            pbc=pbc_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Supercell                                                         #
# ------------------------------------------------------------------ #

def supercell_form(parent, cfg: dict):
    poscar_var = tk.StringVar(value=str(_cfg_get(cfg, "poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var, cfg, "poscar")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    ase_var = _checkbox(parent, "Use ASE backend", _cfg_get(cfg, "by_ase", False))

    def get_args():
        return argparse.Namespace(
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            factors=[int(v.get()) for v in factors_vars],
            by_ase=ase_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Compare                                                           #
# ------------------------------------------------------------------ #

def compare_form(parent, cfg: dict):
    p1_var = tk.StringVar()
    _file_row(parent, "POSCAR 1", p1_var, cfg, "")
    p2_var = tk.StringVar()
    _file_row(parent, "POSCAR 2", p2_var, cfg, "")

    def get_args():
        return argparse.Namespace(poscar1=p1_var.get(), poscar2=p2_var.get())

    return get_args


# ------------------------------------------------------------------ #
#  Merge                                                             #
# ------------------------------------------------------------------ #

def merge_form(parent, cfg: dict):
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)

    tk.Label(parent, text="POSCAR files to merge:").pack(anchor="w", pady=(8, 2))
    listbox = tk.Listbox(parent, height=6)
    listbox.pack(fill=tk.X)

    btn_f = tk.Frame(parent)
    btn_f.pack(fill=tk.X, pady=2)

    def add_file():
        path = filedialog.askopenfilename(
            title="Select POSCAR file",
            filetypes=[("VASP files", "*.vasp"), ("All files", "*.*")],
        )
        if path:
            listbox.insert(tk.END, path)

    def remove_file():
        sel = listbox.curselection()
        if sel:
            listbox.delete(sel[0])

    tk.Button(btn_f, text="+ Add", command=add_file).pack(side=tk.LEFT, padx=2)
    tk.Button(btn_f, text="- Remove", command=remove_file).pack(side=tk.LEFT, padx=2)

    def get_args():
        paths = [listbox.get(i) for i in range(listbox.size())]
        return argparse.Namespace(poscars=paths, outdir=outdir_var.get())

    return get_args


# ------------------------------------------------------------------ #
#  Separate                                                          #
# ------------------------------------------------------------------ #

def separate_form(parent, cfg: dict):
    poscar_var = tk.StringVar()
    _file_row(parent, "POSCAR file", poscar_var, cfg, "")
    outdir_var = tk.StringVar(value=str(_cfg_get(cfg, "outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    _, key_var = _combo(parent, "Group by", ["note", "symbol", "x", "y", "z"], cfg, "separate_key", "note")

    def get_args():
        return argparse.Namespace(
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            key=key_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Registry                                                          #
# ------------------------------------------------------------------ #

FORMS = {
    "Modeling": modeling_form,
    "Import SOFs": import_sofs_form,
    "Count CN": countcn_form,
    "Slice": slice_form,
    "Slice + CN": slice_to_countcn_form,
    "Supercell": supercell_form,
    "Compare": compare_form,
    "Merge": merge_form,
    "Separate": separate_form,
}
