"""Form builders — each returns (frame, get_args_callable)."""

import argparse
import tkinter as tk
from tkinter import ttk, filedialog
from pathlib import Path


# ------------------------------------------------------------------ #
#  Shared widget helpers                                             #
# ------------------------------------------------------------------ #

def _entry(parent, label, cfg, key, default=""):
    """Labelled Entry row — self-packing."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    raw = cfg.get(key, default)
    if isinstance(raw, list):
        raw = " ".join(str(x) for x in raw)
    var = tk.StringVar(value=str(raw))
    tk.Entry(f, textvariable=var).pack(side=tk.LEFT, fill=tk.X, expand=True)
    return var


def _file_row(parent, label, var):
    """Entry + Browse button row."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    tk.Entry(f, textvariable=var).pack(side=tk.LEFT, fill=tk.X, expand=True)
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
        filetypes=[("Config/Data files", "*.toml *.vasp *.csv *.xlsx *.tdb"),
                   ("All files", "*.*")],
    )
    if path:
        var.set(path)


def _browse_dir(var: tk.StringVar):
    path = filedialog.askdirectory(title="Select directory")
    if path:
        var.set(path)


def _checkbox(parent, label, default=False):
    """Create a self-packing checkbox."""
    v = tk.BooleanVar(value=bool(default))
    tk.Checkbutton(parent, text=label, variable=v).pack(anchor="w")
    return v


def _combo(parent, label, values, cfg, key, default=""):
    """Labelled Combobox."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    var = tk.StringVar(value=str(cfg.get(key, default)))
    cb = ttk.Combobox(f, textvariable=var, values=values, state="readonly")
    cb.pack(side=tk.LEFT, fill=tk.X, expand=True)
    return cb, var


def _int_entries(parent, label, count, cfg, key, defaults):
    """Row of N integer entries."""
    f = tk.Frame(parent)
    f.pack(fill=tk.X, pady=2)
    tk.Label(f, text=label, width=16, anchor="e").pack(side=tk.LEFT, padx=(0, 6))
    vals = cfg.get(key, list(defaults))
    vars_ = []
    for i in range(count):
        v = tk.StringVar(value=str(vals[i] if i < len(vals) else defaults[i]))
        tk.Entry(f, textvariable=v, width=5).pack(side=tk.LEFT, padx=1)
        vars_.append(v)
    return vars_


def _parse_temps(raw: str) -> list[float]:
    """Parse temperature string tolerating K suffix, Chinese/ASCII separators."""
    if not raw or not raw.strip():
        return []
    for ch in (",", "，", "；", ";", "、", "\n"):
        raw = raw.replace(ch, " ")
    result = []
    for token in raw.split():
        token = token.strip().rstrip("Kk")
        if token:
            result.append(float(token))
    return result


def _parse_seeds(raw: str) -> list[int] | None:
    """Parse seeds string tolerating various separators. Returns None if empty."""
    if not raw or not raw.strip():
        return None
    for ch in (",", "，", "；", ";", "、", "\n"):
        raw = raw.replace(ch, " ")
    result = []
    for token in raw.split():
        token = token.strip()
        if token:
            result.append(int(token))
    return result or None


# ------------------------------------------------------------------ #
#  SOF editor                                                        #
# ------------------------------------------------------------------ #

def _sof_editor(parent, phase, cfg):
    """Build an editable SOF table for the given phase.

    Returns (editor_frame, collect_fn) where collect_fn() to {site: {elem: fraction}}.
    If phase not found in config, returns (None, lambda: {}).
    """
    phase_cfg = cfg.get(phase.upper(), {})
    if not phase_cfg or "cell" not in phase_cfg:
        return None, lambda: {}

    editor = tk.LabelFrame(parent, text=f"SOFs for {phase.upper()}")
    editor.pack(fill=tk.X, pady=(8, 0))

    rows: list[tuple[str, tk.StringVar, tk.StringVar]] = []
    sum_labels: dict[str, tk.Label] = {}

    def _refresh_sums():
        site_sums: dict[str, float] = {}
        for site, ev, sv in rows:
            try:
                site_sums[site] = site_sums.get(site, 0) + float(sv.get() or 0)
            except ValueError:
                pass
        for site, lbl in sum_labels.items():
            s = site_sums.get(site, 0)
            color = "green" if abs(s - 1.0) < 1e-6 else "red"
            lbl.config(text=f"sum={s:.6f}", fg=color)

    def _add_row(site, before_widget):
        rf = tk.Frame(editor)
        ev = tk.StringVar(value="")
        sv = tk.StringVar(value="0.0")
        e1 = tk.Entry(rf, textvariable=ev, width=8)
        e2 = tk.Entry(rf, textvariable=sv, width=10)
        e1.pack(side=tk.LEFT)
        e2.pack(side=tk.LEFT, padx=4)
        row_ref = [None]  # mutable container for closure
        btn = tk.Button(rf, text="-", width=2,
                        command=lambda: (rf.destroy(), rows.remove(row_ref[0]),
                                         _refresh_sums()))
        btn.pack(side=tk.LEFT)
        rf.pack(fill=tk.X, pady=1, before=before_widget)
        row = (site, ev, sv)
        row_ref[0] = row
        rows.append(row)
        ev.trace_add("write", lambda *a: _refresh_sums())
        sv.trace_add("write", lambda *a: _refresh_sums())

    for site, data in phase_cfg.items():
        if site == "cell":
            continue
        elem, coords = data.get("atoms", ["?", []])
        sofs = data.get("sofs", {})

        header = tk.Label(
            editor, text=f"   {site}   (×{len(coords)}, default: {elem})",
            font=("Arial", 9, "bold"), anchor="w",
        )
        header.pack(fill=tk.X, pady=(6, 2))

        for el, val in sorted(sofs.items()):
            rf = tk.Frame(editor)
            ev = tk.StringVar(value=el)
            sv = tk.StringVar(value=str(val))
            e1 = tk.Entry(rf, textvariable=ev, width=8)
            e2 = tk.Entry(rf, textvariable=sv, width=10)
            e1.pack(side=tk.LEFT)
            e2.pack(side=tk.LEFT, padx=4)
            tk.Button(rf, text="-", width=2,
                      command=lambda r=rf, row=(site, ev, sv):
                          (r.destroy(), rows.remove(row), _refresh_sums())
                      ).pack(side=tk.LEFT)
            rf.pack(fill=tk.X, pady=1)
            row = (site, ev, sv)
            rows.append(row)
            ev.trace_add("write", lambda *a: _refresh_sums())
            sv.trace_add("write", lambda *a: _refresh_sums())

        add_btn = tk.Button(editor, text=f"+ Add to {site}")
        add_btn.configure(command=lambda s=site, b=add_btn: _add_row(s, b))
        add_btn.pack(anchor="w", padx=20)

        sum_lbl = tk.Label(editor, text="sum=0.0000", fg="gray")
        sum_lbl.pack(anchor="e")
        sum_labels[site] = sum_lbl

    _refresh_sums()

    def collect():
        result: dict[str, dict[str, float]] = {}
        for site, ev, sv in rows:
            e = ev.get().strip()
            if not e:
                continue
            try:
                result.setdefault(site, {})[e] = float(sv.get() or 0)
            except ValueError:
                pass
        return result

    return editor, collect


def _link_ase_pbc(ase_var: tk.BooleanVar, pbc_var: tk.BooleanVar):
    """When ASE is checked, auto-check PBC (ASE always uses PBC)."""
    def _sync(*_):
        if ase_var.get():
            pbc_var.set(True)
    ase_var.trace_add("write", _sync)


def _write_sofs_to_config(cfg_path, phase, sofs_by_site):
    """Write edited SOFs into config.toml sections like [PHASE.SITE.sofs]."""
    import re
    phase_upper = phase.upper()
    with open(cfg_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    # Build map: section_header -> {elem: value, ...}
    wanted = {f"{phase_upper}.{site}.sofs": sofs
              for site, sofs in sofs_by_site.items()}

    new_lines = []
    skip_section = False
    seen = set()
    for i, line in enumerate(lines):
        m = re.match(r"^\s*\[([^\]]+)\]\s*$", line)
        if m:
            section = m.group(1).strip()
            skip_section = section in wanted
            if skip_section:
                new_lines.append(line)
                for elem, val in sorted(wanted[section].items()):
                    new_lines.append(f"{elem} = {val:.6g}\n")
                new_lines.append("\n")
                seen.add(section)
                continue
        if skip_section and line.strip() and not line.strip().startswith("#"):
            continue
        new_lines.append(line)

    # Append any section not yet written
    for section, sofs in wanted.items():
        if section not in seen:
            new_lines.append(f"\n[{section}]\n")
            for elem, val in sorted(sofs.items()):
                new_lines.append(f"{elem} = {val:.6g}\n")

    with open(cfg_path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)


# ------------------------------------------------------------------ #
#  Modeling                                                          #
# ------------------------------------------------------------------ #

def modeling_form(parent, cfg: dict):
    import tomllib
    from src.config import normalize_config_keys

    name_var = _entry(parent, "Name", cfg, "name", "modeling")
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    _, phase_var = _combo(parent, "Phase", ["", "FCC", "BCC", "HCP"], cfg, "phase")
    config_var = tk.StringVar(value="config.toml")

    # SOF editor container (rebuilt on phase/config change)
    sof_container = tk.Frame(parent)
    sof_container.pack(fill=tk.X)
    _sof_collect = lambda: {}  # mutable closure

    def _reload_sof(*_):
        nonlocal _sof_collect
        for w in sof_container.winfo_children():
            w.destroy()
        p = phase_var.get().strip()
        cp = config_var.get().strip()
        if not p or not cp or not Path(cp).is_file():
            return
        with open(cp, "rb") as f:
            cfg2 = normalize_config_keys(tomllib.load(f))
        editor, collect = _sof_editor(sof_container, p, cfg2)
        _sof_collect = collect if editor else (lambda: {})

    phase_var.trace_add("write", _reload_sof)
    _file_row(parent, "Config file", config_var)
    config_var.trace_add("write", _reload_sof)

    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    seeds_var = _entry(parent, "Seeds (space-sep)", cfg, "shuffle_seeds", "")
    batch_var = _entry(parent, "Batch size", cfg, "batch_size", "1")
    sqs_var = _checkbox(parent, "Enable SQS", cfg.get("enable_sqs", False))
    iter_var = _entry(parent, "SQS iterations", cfg, "iterations", "10000000")

    def get_args():
        sofs_by_site = _sof_collect()
        phase = phase_var.get()
        cp = config_var.get().strip()
        if sofs_by_site and cp and Path(cp).is_file():
            _write_sofs_to_config(cp, phase, sofs_by_site)
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            phase=phase,
            config=cp or None,
            factors=[int(v.get()) for v in factors_vars],
            seeds=_parse_seeds(seeds_var.get()),
            batch_size=int(batch_var.get() or "1"),
            enable_sqs=sqs_var.get(),
            iterations=int(float(iter_var.get() or "1e7")),
        )

    # Trigger initial load if phase is pre-filled
    if phase_var.get().strip():
        _reload_sof()

    return get_args


# ------------------------------------------------------------------ #
#  Import SOFs                                                       #
# ------------------------------------------------------------------ #

def import_to_modeling_form(parent, cfg: dict):
    csv_var = tk.StringVar()
    _file_row(parent, "CSV/XLSX file", csv_var)
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

    _, output_var = _combo(parent, "Output mode", ["run", "save-config", "print"], cfg, "output", "run")
    outdir_var = tk.StringVar(value="output")
    _dir_row(parent, "Output dir", outdir_var)
    name_var = _entry(parent, "Name", cfg, "name", "modeling")
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    seeds_var = _entry(parent, "Seeds (space-sep)", cfg, "shuffle_seeds", "")
    batch_var = _entry(parent, "Batch size", cfg, "batch_size", "1")
    sqs_var = _checkbox(parent, "Enable SQS", cfg.get("enable_sqs", False))
    iter_var = _entry(parent, "SQS iterations", cfg, "iterations", "10000000")

    def get_args():
        temps = _parse_temps(temps_var.get())
        return argparse.Namespace(
            csv=csv_var.get(),
            config="config.toml",
            phase=phase_var.get(),
            temperatures=temps,
            sublattice_map=map_var.get() or None,
            outdir=outdir_var.get(),
            name=name_var.get(),
            factors=[int(v.get()) for v in factors_vars],
            output=output_var.get(),
            seeds=_parse_seeds(seeds_var.get()),
            batch_size=int(batch_var.get() or "1"),
            enable_sqs=sqs_var.get(),
            iterations=int(float(iter_var.get() or "1e7")),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Count CN                                                          #
# ------------------------------------------------------------------ #

def countcn_form(parent, cfg: dict):
    name_var = _entry(parent, "Name", cfg, "name", "countcn")
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    mult_var = _entry(parent, "Cutoff multiplier", cfg, "cutoff_mult", "1.1")
    par_var = _entry(parent, "Parallel workers", cfg, "parallel", "2")
    ase_var = _checkbox(parent, "Use ASE backend", cfg.get("by_ase", False))
    pbc_var = _checkbox(parent, "Periodic boundary (PBC)", cfg.get("pbc", False))
    _link_ase_pbc(ase_var, pbc_var)

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
    name_var = _entry(parent, "Name", cfg, "name", "slice")
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
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
#  Slice to CountCN                                                        #
# ------------------------------------------------------------------ #

def slice_to_countcn_form(parent, cfg: dict):
    name_var = _entry(parent, "Name", cfg, "name", "slice-to-countcn")
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    miller_vars = _int_entries(parent, "Miller index", 3, cfg,
                               "slice_direction", (0, 0, 1))
    ase_var = _checkbox(parent, "Use ASE backend", cfg.get("by_ase", False))
    pbc_var = _checkbox(parent, "Periodic boundary (PBC)", cfg.get("pbc", False))
    _link_ase_pbc(ase_var, pbc_var)

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            miller_index=[int(v.get()) for v in miller_vars],
            pbc=pbc_var.get(),
            by_ase=ase_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Supercell                                                         #
# ------------------------------------------------------------------ #

def supercell_form(parent, cfg: dict):
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    ase_var = _checkbox(parent, "Use ASE backend", cfg.get("by_ase", False))

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
    _file_row(parent, "POSCAR 1", p1_var)
    p2_var = tk.StringVar()
    _file_row(parent, "POSCAR 2", p2_var)

    def get_args():
        return argparse.Namespace(poscar1=p1_var.get(), poscar2=p2_var.get())

    return get_args


# ------------------------------------------------------------------ #
#  Merge                                                             #
# ------------------------------------------------------------------ #

def merge_form(parent, cfg: dict):
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
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
#  Thermo                                                            #
# ------------------------------------------------------------------ #

def thermo_form(parent, cfg: dict):
    name_var = _entry(parent, "Name", cfg, "name", "thermo")
    data_var = tk.StringVar()
    _file_row(parent, "Data file (xlsx/csv)", data_var)
    tdb_var = tk.StringVar()
    _file_row(parent, "TDB file", tdb_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)

    def get_args():
        return argparse.Namespace(
            name=name_var.get(),
            data=data_var.get(),
            tdb=tdb_var.get(),
            outdir=outdir_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Separate                                                          #
# ------------------------------------------------------------------ #

def separate_form(parent, cfg: dict):
    poscar_var = tk.StringVar()
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
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
#  Config editor                                                     #
# ------------------------------------------------------------------ #

def config_form(parent, cfg: dict):
    """Edit all config.toml parameters directly."""
    import tomllib
    from src.config import normalize_config_keys

    config_var = tk.StringVar(value="config.toml")
    _file_row(parent, "Config file", config_var)

    tk.Label(parent, text="General", font=("Arial", 10, "bold"),
             anchor="w").pack(fill=tk.X, pady=(10, 2))
    name_var = _entry(parent, "Name", cfg, "name", "modeling")
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    _, phase_var = _combo(parent, "Phase", ["", "FCC", "BCC", "HCP"], cfg, "phase")

    # SOF editor container (rebuilt on phase/config change)
    sof_container = tk.Frame(parent)
    sof_container.pack(fill=tk.X)
    _sof_collect = lambda: {}

    def _reload_sof(*_):
        nonlocal _sof_collect
        for w in sof_container.winfo_children():
            w.destroy()
        p = phase_var.get().strip()
        cp = config_var.get().strip()
        if not p or not cp or not Path(cp).is_file():
            _sof_collect = lambda: {}
            return
        with open(cp, "rb") as f:
            cfg2 = normalize_config_keys(tomllib.load(f))
        editor, collect = _sof_editor(sof_container, p, cfg2)
        _sof_collect = collect if editor else (lambda: {})

    phase_var.trace_add("write", _reload_sof)
    config_var.trace_add("write", _reload_sof)

    tk.Label(parent, text="Supercell & Modeling", font=("Arial", 10, "bold"),
             anchor="w").pack(fill=tk.X, pady=(10, 2))
    factors_vars = _int_entries(parent, "Supercell factors", 3, cfg,
                                "supercell_factors", (3, 3, 3))
    seeds_var = _entry(parent, "Seeds (space-sep)", cfg, "shuffle_seeds", "")
    batch_var = _entry(parent, "Batch size", cfg, "batch_size", "1")
    sqs_var = _checkbox(parent, "Enable SQS", cfg.get("enable_sqs", False))
    iter_var = _entry(parent, "SQS iterations", cfg, "iterations", "10000000")

    tk.Label(parent, text="Analysis", font=("Arial", 10, "bold"),
             anchor="w").pack(fill=tk.X, pady=(10, 2))
    mult_var = _entry(parent, "Cutoff multiplier", cfg, "cutoff_mult", "1.1")
    par_var = _entry(parent, "Parallel workers", cfg, "parallel", "2")
    ase_var = _checkbox(parent, "Use ASE backend", cfg.get("by_ase", False))
    pbc_var = _checkbox(parent, "Periodic boundary (PBC)", cfg.get("pbc", False))
    miller_vars = _int_entries(parent, "Miller index", 3, cfg,
                               "slice_direction", (0, 0, 1))

    # Trigger initial load if phase is pre-filled
    if phase_var.get().strip():
        _reload_sof()

    def get_args():
        sofs_by_site = _sof_collect()
        cp = config_var.get().strip()
        phase = phase_var.get()
        if sofs_by_site and cp and Path(cp).is_file():
            _write_sofs_to_config(cp, phase, sofs_by_site)
        return argparse.Namespace(
            config_path=cp,
            name=name_var.get(),
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            phase=phase,
            factors=[int(v.get()) for v in factors_vars],
            seeds=_parse_seeds(seeds_var.get()),
            batch_size=int(batch_var.get() or "1"),
            enable_sqs=sqs_var.get(),
            iterations=int(float(iter_var.get() or "1e7")),
            cutoff_mult=float(mult_var.get() or "1.1"),
            parallel=int(par_var.get() or "2"),
            by_ase=ase_var.get(),
            pbc=pbc_var.get(),
            slice_direction=[int(v.get()) for v in miller_vars],
        )

    return get_args


# ------------------------------------------------------------------ #
#  Surface                                                           #
# ------------------------------------------------------------------ #

def surface_form(parent, cfg: dict):
    poscar_var = tk.StringVar(value=str(cfg.get("poscar", "")))
    _file_row(parent, "POSCAR file", poscar_var)
    outdir_var = tk.StringVar(value=str(cfg.get("outdir", "output")))
    _dir_row(parent, "Output dir", outdir_var)
    miller_vars = _int_entries(parent, "Miller index", 3, cfg,
                               "slice_direction", (0, 0, 1))
    layers_var = _entry(parent, "Layers (N)", cfg, "layers", "3")
    vacuum_var = _entry(parent, "Vacuum (Å)", cfg, "vacuum", "15.0")
    fix_layers_var = _entry(parent, "Fix layers (blank=auto)", cfg, "fix_layers", "")
    fix_z_var = _checkbox(parent, "Fix z-only (TTF)", cfg.get("fix_z_only", False))

    def get_args():
        fix_raw = fix_layers_var.get().strip()
        return argparse.Namespace(
            poscar=poscar_var.get(),
            outdir=outdir_var.get(),
            miller=[int(v.get()) for v in miller_vars],
            layers=int(layers_var.get() or "3"),
            vacuum=float(vacuum_var.get() or "15.0"),
            fix_layers=int(fix_raw) if fix_raw else None,
            fix_z_only=fix_z_var.get(),
        )

    return get_args


# ------------------------------------------------------------------ #
#  Registry                                                          #
# ------------------------------------------------------------------ #

DESCRIPTIONS: dict[str, str] = {
    "Config": "View and edit all config.toml parameters. Load/Save configuration.",
    "Modeling": "Generate supercell and allocate atoms based on SOFs.",
    "Import to Model": (
        "Import site-fraction data from ThermoCalc/Pandat CSV/XLSX\n"
        "and run modeling at selected temperatures."
    ),
    "Count CN": (
        "Count coordination numbers and nearest-neighbor pairs.\n"
        "Supports KDTree and ASE backends, optional PBC."
    ),
    "Slice": "Slice a structure normal to a Miller index and export layers.",
    "Slice to CountCN": "Slice a structure, then count CN for each layer separately.",
    "Supercell": (
        "Expand a unit cell along basis vectors.\n"
        "Optionally use ASE's make_supercell."
    ),
    "Compare": "Compare two POSCAR files (cell, symbols, atom positions).",
    "Merge": "Combine multiple POSCAR files into one.",
    "Separate": "Split a POSCAR file by element, coordinate, or note.",
    "Thermo": (
        "Calculate configurational entropy (Sconf) and Gibbs free energy\n"
        "(DeltaG) from SOF data and TDB thermodynamic database."
    ),
    "About": "Developer information, license, and project links.",
    "Surface": (
        "Generate asymmetric surface slabs from a bulk POSCAR.\n"
        "Select Miller index, number of layers (N and N+1 output),\n"
        "vacuum thickness, and optional bottom-layer fixing."
    ),
}

FORMS = {
    "Config": config_form,
    "Modeling": modeling_form,
    "Import to Model": import_to_modeling_form,
    "Count CN": countcn_form,
    "Slice": slice_form,
    "Slice to CountCN": slice_to_countcn_form,
    "Surface": surface_form,
    "Supercell": supercell_form,
    "Compare": compare_form,
    "Merge": merge_form,
    "Separate": separate_form,
    "Thermo": thermo_form,
}
