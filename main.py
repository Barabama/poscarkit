import sys
from src.cli.poscarkit import main as main_cli
from src.cli.poscarkit_interact import main as main_interact

if __name__ == "__main__":
    # Check if there are any arguments
    if len(sys.argv) > 1:
        # If there are arguments, use the CLI mode
        main_cli()
    else:
        # If there are no arguments, use the interactive mode
        main_interact()

# 编译指令
# nuitka --standalone --onefile --output-dir=dist --jobs=4 --lto=yes `
# --enable-plugin=tk-inter --enable-plugin=no-qt --windows-console-mode=disable `
# --windows-icon-from-ico="icon.ico" `
# --onefile-no-compression `
# --enable-plugin=upx --upx-binary="D:\\Programs\\upx-5.0.2-win64\\upx.exe" `
# --output-filename=poscarkit-0.9.0.exe `
# --file-version=0.9.0 `
# --copyright="(C) 2025 MCMF, Fuzhou University" `
# main.py
