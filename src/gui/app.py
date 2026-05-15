"""Tkinter GUI — sidebar navigation + dynamic form + log area."""

import logging
import sys
import threading
import tkinter as tk
from tkinter import ttk

from src.gui.forms import FORMS
from src.config import VERSION

FUNC_NAMES = [
    ("Modeling", "#4A90D9"),
    ("Import SOFs", "#50B86C"),
    ("Count CN", "#E8A838"),
    ("Slice", "#D94A4A"),
    ("Slice + CN", "#9B59B6"),
    ("Supercell", "#E67E22"),
    ("Compare", "#1ABC9C"),
    ("Merge", "#3498DB"),
    ("Separate", "#95A5A6"),
]


class _LogHandler(logging.Handler):
    """Redirect logging to a tk.Text widget (thread-safe)."""

    def __init__(self, widget: tk.Text):
        super().__init__()
        self.widget = widget
        self.setFormatter(logging.Formatter("%(asctime)s[%(levelname)s]%(message)s",
                                             datefmt="%Y-%m-%d %H:%M:%S"))

    def emit(self, record):
        msg = self.format(record) + "\n"
        self.widget.after(0, self._append, msg)

    def _append(self, msg):
        self.widget.insert(tk.END, msg)
        self.widget.see(tk.END)


class PoscaKitGUI:
    """Main GUI application window."""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title(f"POSCARKIT {VERSION.split(':')[1].split('|')[0].strip()}")
        self.root.geometry("960x780")
        self.root.minsize(800, 500)

        # Styling
        self.style = ttk.Style()
        self.style.theme_use("clam")

        # Config
        self._cfg: dict = {}
        self._load_config()

        # State (init before building widgets that reference these)
        self._current_form = None
        self._func_buttons: dict[str, tk.Button] = {}
        self._running = False

        # Layout frames
        self._build_sidebar()
        self._build_main_area()

        # Log redirect
        self._log_handler = _LogHandler(self._log_area)
        logging.getLogger().addHandler(self._log_handler)
        logging.getLogger().setLevel(logging.INFO)
        logging.info("GUI started. Select a function from the sidebar.")

        # Select first function by default
        self._switch_form(FUNC_NAMES[0][0])

    # ------------------------------------------------------------------ #
    #  Sidebar                                                           #
    # ------------------------------------------------------------------ #

    def _build_sidebar(self):
        sidebar = tk.Frame(self.root, bg="#2C3E50", width=150)
        sidebar.pack(side=tk.LEFT, fill=tk.Y)
        sidebar.pack_propagate(False)

        tk.Label(
            sidebar, text="POSCARKIT", bg="#2C3E50", fg="white",
            font=("Arial", 12, "bold"), pady=8,
        ).pack()

        for name, color in FUNC_NAMES:
            btn = tk.Button(
                sidebar, text=name, anchor="w", relief=tk.FLAT,
                bg="#34495E", fg="white", font=("Arial", 10),
                activebackground=color, activeforeground="white",
                padx=12, pady=6, borderwidth=0,
                command=lambda n=name: self._switch_form(n),
            )
            btn.pack(fill=tk.X)
            self._func_buttons[name] = btn

    # ------------------------------------------------------------------ #
    #  Main area (form + log)                                            #
    # ------------------------------------------------------------------ #

    def _build_main_area(self):
        right = tk.Frame(self.root)
        right.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

        # PanedWindow: form on top, log on bottom (user can drag to resize)
        pane = tk.PanedWindow(right, orient=tk.VERTICAL, sashrelief=tk.RAISED, sashwidth=4)
        pane.pack(fill=tk.BOTH, expand=True)

        # Top: form area inside a canvas + scrollbar
        form_outer = tk.Frame(pane)
        self._title_label = tk.Label(
            form_outer, text="", font=("Arial", 14, "bold"), anchor="w", pady=6, padx=12,
        )
        self._title_label.pack(fill=tk.X)

        form_canvas = tk.Canvas(form_outer, highlightthickness=0)
        form_scroll = ttk.Scrollbar(form_outer, orient="vertical", command=form_canvas.yview)
        self._form_frame = tk.Frame(form_canvas, padx=12, pady=6)
        _fw = form_canvas.create_window((0, 0), window=self._form_frame, anchor="nw")
        self._form_frame.bind(
            "<Configure>",
            lambda e: form_canvas.configure(scrollregion=form_canvas.bbox("all")),
        )
        form_canvas.bind("<Configure>", lambda e, w=_fw: form_canvas.itemconfig(w, width=e.width))
        form_canvas.configure(yscrollcommand=form_scroll.set)
        form_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        form_scroll.pack(side=tk.RIGHT, fill=tk.Y)

        pane.add(form_outer, stretch="always")

        # Bottom: log area
        log_frame = tk.LabelFrame(pane, text="Log", padx=4, pady=2)
        pane.add(log_frame, stretch="always")

        self._log_area = tk.Text(
            log_frame, height=12, wrap=tk.WORD, state=tk.NORMAL,
            font=("Consolas", 9), bg="#1E1E1E", fg="#D4D4D4",
        )
        scroll = ttk.Scrollbar(log_frame, command=self._log_area.yview)
        self._log_area.configure(yscrollcommand=scroll.set)
        self._log_area.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)

    # ------------------------------------------------------------------ #
    #  Form switching                                                    #
    # ------------------------------------------------------------------ #

    def _switch_form(self, name: str):
        if self._running:
            logging.warning("A task is running. Please wait for it to finish.")
            return

        # Highlight active button
        for n, btn in self._func_buttons.items():
            btn.configure(bg="#34495E" if n != name else _color_for_name(name))

        self._title_label.configure(text=name)

        # Clear old form
        for w in self._form_frame.winfo_children():
            w.destroy()

        # Build new form
        builder = FORMS.get(name)
        if builder is None:
            tk.Label(self._form_frame, text="(not implemented yet)").pack()
            return

        get_args = builder(self._form_frame, self._cfg)

        # Run button
        btn_frame = tk.Frame(self._form_frame)
        btn_frame.pack(fill=tk.X, pady=(12, 0))

        self._run_btn = tk.Button(
            btn_frame, text=f"Run {name}", bg="#27AE60", fg="white",
            font=("Arial", 11, "bold"), padx=24, pady=4,
            command=lambda: self._run_task(name, get_args),
        )
        self._run_btn.pack(side=tk.LEFT, padx=(0, 8))

        # Save checkbox + Load button
        self._save_var = tk.BooleanVar(value=False)
        tk.Checkbutton(
            btn_frame, text="Save as defaults", variable=self._save_var,
        ).pack(side=tk.LEFT, padx=8)

        tk.Button(
            btn_frame, text="Load from config",
            command=lambda: self._load_config_and_rebuild(name),
        ).pack(side=tk.LEFT, padx=8)

    def _load_config_and_rebuild(self, name: str):
        self._load_config()
        self._switch_form(name)

    # ------------------------------------------------------------------ #
    #  Run task in background thread                                     #
    # ------------------------------------------------------------------ #

    def _run_task(self, name: str, get_args):
        if self._running:
            return

        try:
            args = get_args()
        except Exception as e:
            logging.error(f"Invalid parameters: {e}")
            return

        self._running = True
        self._run_btn.configure(text="Running...", state=tk.DISABLED, bg="#95A5A6")

        def target():
            try:
                # Lazy-import the handler
                from src.cli.poscarkit import (
                    cmd_modeling, cmd_countcn, cmd_slice,
                    cmd_slice_to_countcn, cmd_supercell,
                    cmd_compare, cmd_merge, cmd_separate,
                )
                handlers = {
                    "Modeling": cmd_modeling,
                    "Count CN": cmd_countcn,
                    "Slice": cmd_slice,
                    "Slice + CN": cmd_slice_to_countcn,
                    "Supercell": cmd_supercell,
                    "Compare": cmd_compare,
                    "Merge": cmd_merge,
                    "Separate": cmd_separate,
                }
                handler = handlers.get(name)
                if handler:
                    code = handler(args)
                    if code == 0 and self._save_var.get():
                        self._save_to_config(args)
                    logging.info(f"Task '{name}' finished (exit {code}).")
                else:
                    logging.error(f"No handler for '{name}'")
            except Exception:
                logging.exception(f"Task '{name}' failed")
            finally:
                self.root.after(0, self._task_done, name)

        threading.Thread(target=target, daemon=True).start()

    def _task_done(self, name: str):
        self._running = False
        if hasattr(self, "_run_btn"):
            self._run_btn.configure(
                text=f"Run {name}", state=tk.NORMAL, bg="#27AE60"
            )

    # ------------------------------------------------------------------ #
    #  Config I/O                                                        #
    # ------------------------------------------------------------------ #

    def _load_config(self):
        import tomllib
        from pathlib import Path
        from src.config import normalize_config_keys

        cfg_path = Path("config.toml")
        if cfg_path.is_file():
            with open(cfg_path, "rb") as f:
                self._cfg = normalize_config_keys(tomllib.load(f))
        else:
            self._cfg = {}

    def _save_to_config(self, args):
        """Naive key=value update of config.toml from namespace fields."""
        from pathlib import Path
        from src.config import normalize_config_keys

        # Map arg names to config keys
        key_map = {
            "name": "name", "poscar": "poscar", "outdir": "outdir",
            "phase": "phase", "factors": "supercell_factors",
            "seeds": "shuffle_seeds", "batch_size": "batch_size",
            "enable_sqs": "enable_sqs", "cutoff_mult": "cutoff_mult",
            "by_ase": "by_ase", "pbc": "pbc",
        }

        updates = {}
        for arg_k, cfg_k in key_map.items():
            val = getattr(args, arg_k, None)
            if val is not None:
                if isinstance(val, tuple):
                    val = list(val)
                updates[cfg_k] = val

        # Read current file, replace matching lines
        cfg_path = Path("config.toml")
        if not cfg_path.is_file():
            # Write a minimal config
            with open(cfg_path, "w", encoding="utf-8") as f:
                for k, v in updates.items():
                    f.write(f"{k} = {v!r}\n")
            self._load_config()
            logging.info("Config saved.")
            return

        import re
        with open(cfg_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        written = set()
        new_lines = []
        for line in lines:
            m = re.match(r"^(\w+)\s*=", line.strip())
            if m and m.group(1) in updates:
                key = m.group(1)
                val = updates[key]
                new_lines.append(f"{key} = {_toml_value(val)}\n")
                written.add(key)
                continue
            new_lines.append(line)

        # Append keys not yet written
        for k, v in updates.items():
            if k not in written:
                new_lines.append(f"{k} = {_toml_value(v)}\n")

        with open(cfg_path, "w", encoding="utf-8") as f:
            f.writelines(new_lines)

        self._load_config()
        logging.info("Config saved.")

    # ------------------------------------------------------------------ #
    #  Public API                                                        #
    # ------------------------------------------------------------------ #

    def run(self):
        """Enter the Tk main loop."""
        self.root.mainloop()


# -- helpers ----------------------------------------------------------- #


def _color_for_name(name: str) -> str:
    for n, c in FUNC_NAMES:
        if n == name:
            return c
    return "#34495E"


def _toml_value(val) -> str:
    if isinstance(val, bool):
        return "true" if val else "false"
    if isinstance(val, list):
        return "[" + ", ".join(str(v) for v in val) + "]"
    if isinstance(val, str):
        return f'"{val}"'
    return str(val)
