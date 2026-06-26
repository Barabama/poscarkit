# main.py

import logging
import os
import sys
import tomllib
import traceback

from GTCalc import GTCalculator


VERSION = "0.1.0"
INFO_EXEC = f"""
------------------------------------------
G(T)-Calculator (ver {VERSION})
------------------------------------------
功能描述:
  给定某相的亚晶格的占位分数(有序), 输出该相的G-T数据(有序/无序).

使用说明:
  1. 配置 config.toml 文件, 指定相名称以及亚晶格的占位分数.
  2. 运行程序.
  3. 提供 .tdb 文件, 输出 G-T 数据.
------------------------------------------
"""
DEF_CFG = """
# Built in default config.toml

# tdb file path
tdb = "BCC.tdb"
# phase_name
phase = "BCC"

# site of fractions for phase
[BCC.1a.sofs]
Fe = 0.5
Ni = 0.5
[BCC.1b.sofs]
Fe = 0.5
Ni = 0.5
"""

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")


def read_config():
    cfg_path = os.path.join(os.path.dirname(sys.argv[0]), "config.toml")
    if not os.path.isfile(cfg_path):
        logging.error("config.toml missing, building default config.toml")
        with open(cfg_path, "w", encoding="utf-8") as f:
            f.write(DEF_CFG)

    with open(cfg_path, "rb") as f:
        cfg = tomllib.load(f)
    return cfg


def handle_tdb(cfg: dict, tdb: str = ""):
    tdb_cfg = cfg.get("tdb", "")
    while True:
        tdb = tdb or input("Enter tdb file path(default to config.tdb): ") or tdb_cfg
        if not os.path.isfile(tdb):
            logging.warning(f"TDB {tdb} not found, please try again.")
            tdb_cfg = ""
            tdb = ""

        return tdb


def handle_phase(cfg: dict, phase: str = ""):
    phase_cfg = cfg.get("phase", "")
    while True:
        phase = phase or input("Enter phase name(default to config.phase): ").upper() or phase_cfg
        if phase not in cfg:
            logging.warning(f"Phase {phase} not found in config, please try again.")
            phase = ""
        try:
            return phase, {s.upper(): d["sofs"] for s, d in cfg[phase].items()}
        except KeyError:
            logging.error(f"Phase {phase} missing sofs, please try again.")


def main():
    print(INFO_EXEC)
    tdb = sys.argv[1] if len(sys.argv) > 1 else ""
    phase = sys.argv[2] if len(sys.argv) > 2 else ""
    while True:
        try:
            cfg = read_config()
            tdb = handle_tdb(cfg, tdb)
            phase, site_fracts = handle_phase(cfg, phase)
            gtc = GTCalculator(tdb)
            gtc.handle(site_fracts, phase)
            logging.info(f"G-T data for {phase} saved to /output/ folder.")

            # Reset
            logging.info("Resetting and re-reading config file...")
            tdb = ""
            phase = ""
        except KeyboardInterrupt:
            print("\n")
            logging.info("Exiting...")
            sys.exit(0)
        except Exception as e:
            logging.critical(e)
            traceback.print_exc()
            input("Press Enter to exit...")
            sys.exit(1)


if __name__ == "__main__":
    main()
