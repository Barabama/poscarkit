import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        from src.cli.poscarkit import main as main_cli
        sys.exit(main_cli())
    else:
        from src.gui.app import PoscaKitGUI
        PoscaKitGUI().run()
