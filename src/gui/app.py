"""Tkinter GUI — sidebar navigation + dynamic form + log area."""

import logging
import threading
import tkinter as tk
from collections.abc import Callable
from tkinter import ttk

from src.gui.forms import FORMS, DESCRIPTIONS
from src.config import VERSION

FUNC_NAMES = [
    ("Config", "#7F8C8D"),
    ("Modeling", "#4A90D9"),
    ("Import to Model", "#50B86C"),
    ("Count CN", "#E8A838"),
    ("Slice", "#D94A4A"),
    ("Slice to CountCN", "#9B59B6"),
    ("Thermo", "#1ABC9C"),
    ("Supercell", "#E67E22"),
    ("Compare", "#8E44AD"),
    ("Merge", "#3498DB"),
    ("Separate", "#95A5A6"),
]


class _LogHandler(logging.Handler):
    """Redirect logging to a tk.Text widget (thread-safe)."""

    _MAX_LOG_LINES = 5000

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
        line_count = int(self.widget.index("end-1c").split(".")[0])
        if line_count > self._MAX_LOG_LINES:
            self.widget.delete("1.0", f"{line_count - self._MAX_LOG_LINES}.0")
        self.widget.see(tk.END)


class PoscaKitGUI:
    """Main GUI application window."""

    def __init__(self):
        self.root = tk.Tk()
        self.root.title(f"POSCARKIT {VERSION.split(':')[1].split('|')[0].strip()}")
        self.root.geometry("960x820")
        self.root.minsize(800, 500)

        # Styling
        self.style = ttk.Style()
        self.style.theme_use("clam")

        # Config
        self._cfg: dict = {}
        self._cfg_path: str = "config.toml"
        self._load_config()

        # State (init before building widgets that reference these)
        self._current_form = None
        self._func_buttons: dict[str, tk.Button] = {}
        self._running = False
        self._form_cache: dict[str, tuple[tk.Frame, Callable]] = {}
        self._form_state: dict[str, dict] = {}

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

        # PanedWindow: form on top, log on bottom (draggable divider)
        pane = tk.PanedWindow(right, orient=tk.VERTICAL, sashrelief=tk.RAISED, sashwidth=4)
        pane.pack(fill=tk.BOTH, expand=True)

        # Top: form area with scrollbar
        form_outer = tk.Frame(pane)
        self._title_label = tk.Label(
            form_outer, text="", font=("Arial", 14, "bold"), anchor="w", pady=6, padx=12,
        )
        self._title_label.pack(fill=tk.X)
        self._desc_label = tk.Label(
            form_outer, text="", font=("Arial", 9), fg="gray",
            anchor="w", justify=tk.LEFT, padx=12,
        )
        self._desc_label.pack(fill=tk.X)

        canvas = tk.Canvas(form_outer, highlightthickness=0)
        scrollbar = ttk.Scrollbar(form_outer, orient="vertical", command=canvas.yview)
        self._form_frame = tk.Frame(canvas, padx=12, pady=6)
        self._canvas = canvas
        self._win_id = canvas.create_window((0, 0), window=self._form_frame, anchor="nw")

        def _on_form_configure(event):
            bbox = canvas.bbox("all")
            if bbox is None:
                return
            canvas.configure(scrollregion=bbox)
            if bbox[3] <= canvas.winfo_height():
                scrollbar.pack_forget()
            else:
                scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self._form_frame.bind("<Configure>", _on_form_configure)

        def _on_canvas_configure(event):
            canvas.itemconfig(self._win_id, width=event.width)
            # Re-check scrollbar visibility
            if canvas.bbox("all") and canvas.bbox("all")[3] <= event.height:
                scrollbar.pack_forget()
            else:
                scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        canvas.bind("<Configure>", _on_canvas_configure)

        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Mousewheel scrolling — bind_all so child widgets don't swallow the event
        def _on_mousewheel(event):
            w = event.widget
            try:
                while w is not None:
                    if w is canvas:
                        break
                    w = w.master
                else:
                    return
            except (AttributeError, tk.TclError):
                return
            bbox = canvas.bbox("all")
            if bbox and bbox[3] > canvas.winfo_height():
                canvas.yview_scroll(int(-1 * (event.delta / 120)), "units")

        self.root.bind_all("<MouseWheel>", _on_mousewheel)

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
        self._desc_label.configure(text=DESCRIPTIONS.get(name, ""))

        # Hide current form
        if self._current_form and self._current_form in self._form_cache:
            frame, _ = self._form_cache[self._current_form]
            frame.pack_forget()

        self._current_form = name

        # Show cached form or build new one
        if name in self._form_cache:
            frame, get_args = self._form_cache[name]
            frame.pack(fill=tk.BOTH, expand=True)
        else:
            frame = tk.Frame(self._form_frame, padx=12, pady=6)
            builder = FORMS.get(name)
            if builder is None:
                tk.Label(frame, text="(not implemented yet)").pack()
                frame.pack(fill=tk.BOTH, expand=True)
                self._form_cache[name] = (frame, lambda: None)
                return

            get_args = builder(frame, self._cfg)

            # Action buttons
            btn_frame = tk.Frame(frame)
            btn_frame.pack(fill=tk.X, pady=(12, 0))

            if name == "Config":
                run_btn = tk.Button(
                    btn_frame, text="Save Config", bg="#2980B9", fg="white",
                    font=("Arial", 11, "bold"), padx=24, pady=4,
                    command=lambda: self._save_config_form(get_args),
                )
                run_btn.pack(side=tk.LEFT, padx=(0, 8))
            else:
                run_btn = tk.Button(
                    btn_frame, text=f"Run {name}", bg="#27AE60", fg="white",
                    font=("Arial", 11, "bold"), padx=24, pady=4,
                    command=lambda: self._run_task(name, get_args),
                )
                run_btn.pack(side=tk.LEFT, padx=(0, 8))

            tk.Button(
                btn_frame, text="Load from config",
                command=lambda n=name: self._reload_form(n),
            ).pack(side=tk.LEFT, padx=8)

            frame.pack(fill=tk.BOTH, expand=True)
            self._form_cache[name] = (frame, get_args)
            self._form_state[name] = {"run_btn": run_btn}

        # Update references for _run_task / _task_done
        state = self._form_state.get(name, {})
        self._run_btn = state.get("run_btn")

        # Reset scroll to top when switching forms
        self._canvas.yview_moveto(0)

    def _reload_form(self, name: str):
        """Reload config and rebuild a specific form."""
        if name in self._form_cache:
            _, get_args = self._form_cache[name]
            try:
                args = get_args()
                cp = getattr(args, "config_path", None)
                if cp:
                    self._cfg_path = cp
            except Exception:
                pass
            frame, _ = self._form_cache[name]
            frame.destroy()
            del self._form_cache[name]
        if name in self._form_state:
            del self._form_state[name]
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
        if self._run_btn:
            self._run_btn.configure(text="Running...", state=tk.DISABLED, bg="#95A5A6")

        def target():
            try:
                from src.cli.poscarkit import (
                    cmd_modeling, cmd_countcn, cmd_slice,
                    cmd_slice_to_countcn, cmd_supercell,
                    cmd_compare, cmd_merge, cmd_separate,
                    cmd_import_to_model, cmd_thermo,
                )
                handlers = {
                    "Modeling": cmd_modeling,
                    "Import to Model": cmd_import_to_model,
                    "Count CN": cmd_countcn,
                    "Slice": cmd_slice,
                    "Slice to CountCN": cmd_slice_to_countcn,
                    "Thermo": cmd_thermo,
                    "Supercell": cmd_supercell,
                    "Compare": cmd_compare,
                    "Merge": cmd_merge,
                    "Separate": cmd_separate,
                }
                handler = handlers.get(name)
                if handler:
                    code = handler(args)
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
        if self._run_btn is not None:
            self._run_btn.configure(
                text=f"Run {name}", state=tk.NORMAL, bg="#27AE60"
            )

    def _save_config_form(self, get_args):
        """Save all Config form fields to the config file."""
        try:
            args = get_args()
        except Exception as e:
            logging.error(f"Invalid parameters: {e}")
            return
        config_path = getattr(args, "config_path", None) or "config.toml"
        self._cfg_path = config_path
        self._save_to_config(args)
        logging.info("Config saved to: " + config_path)

    # ------------------------------------------------------------------ #
    #  Config I/O                                                        #
    # ------------------------------------------------------------------ #

    def _load_config(self):
        import tomllib
        from pathlib import Path
        from src.config import normalize_config_keys, DEFAULT_CONFIG

        cfg_path = Path(self._cfg_path)
        if not cfg_path.is_file():
            with open(cfg_path, "w", encoding="utf-8") as f:
                f.write(DEFAULT_CONFIG)
            logging.info(f"Generated default {cfg_path}")

        try:
            with open(cfg_path, "rb") as f:
                self._cfg = normalize_config_keys(tomllib.load(f))
        except Exception as e:
            logging.error(f"Failed to parse {cfg_path}: {e}")
            self._cfg = {}

    def _save_to_config(self, args):
        """Update top-level keys in config file from namespace fields."""
        from pathlib import Path
        import re

        key_map = {
            "name": "name", "poscar": "poscar", "outdir": "outdir",
            "phase": "phase", "factors": "supercell_factors",
            "seeds": "shuffle_seeds", "batch_size": "batch_size",
            "enable_sqs": "enable_sqs", "iterations": "iterations",
            "cutoff_mult": "cutoff_mult", "parallel": "parallel",
            "by_ase": "by_ase", "pbc": "pbc",
            "slice_direction": "slice_direction",
        }

        updates = {}
        for arg_k, cfg_k in key_map.items():
            val = getattr(args, arg_k, None)
            if val is None:
                continue
            if isinstance(val, tuple):
                val = list(val)
            updates[cfg_k] = val

        cfg_path = Path(self._cfg_path)
        if not cfg_path.is_file():
            with open(cfg_path, "w", encoding="utf-8") as f:
                for k, v in updates.items():
                    f.write(f"{k} = {_toml_value(v)}\n")
            self._load_config()
            logging.info("Config saved.")
            return

        with open(cfg_path, "r", encoding="utf-8") as f:
            lines = f.readlines()

        written = set()
        new_lines = []
        first_section_idx = None

        for line in lines:
            stripped = line.strip()
            # Track where the first [section] header is
            if first_section_idx is None and re.match(r"^\[", stripped):
                first_section_idx = len(new_lines)

            # Match uncommented top-level key = value (only before first section or if already top-level)
            m = re.match(r"^(\w+)\s*=", stripped)
            if m and m.group(1) in updates and (first_section_idx is None):
                key = m.group(1)
                new_lines.append(f"{key} = {_toml_value(updates[key])}\n")
                written.add(key)
                continue

            # Also match commented-out keys and uncomment them with new value
            cm = re.match(r"^#\s*(\w+)\s*=", stripped)
            if cm and cm.group(1) in updates and cm.group(1) not in written and (first_section_idx is None):
                key = cm.group(1)
                new_lines.append(f"{key} = {_toml_value(updates[key])}\n")
                written.add(key)
                continue

            new_lines.append(line)

        # Insert remaining keys before the first section header
        remaining = [f"{k} = {_toml_value(v)}\n" for k, v in updates.items() if k not in written]
        if remaining:
            insert_at = first_section_idx if first_section_idx is not None else len(new_lines)
            for line in remaining:
                new_lines.insert(insert_at, line)
                insert_at += 1
            if first_section_idx is not None:
                new_lines.insert(insert_at, "\n")

        with open(cfg_path, "w", encoding="utf-8") as f:
            f.writelines(new_lines)

        self._load_config()
        logging.info(f"Config saved. Updated keys: {list(updates.keys())}")

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
        parts = []
        for v in val:
            if isinstance(v, bool):
                parts.append("true" if v else "false")
            elif isinstance(v, str):
                parts.append(f'"{v}"')
            else:
                parts.append(str(v))
        return "[" + ", ".join(parts) + "]"
    if isinstance(val, str):
        return f'"{val}"'
    return str(val)
