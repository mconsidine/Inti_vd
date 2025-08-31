# -*- mode: python ; coding: utf-8 -*-
block_cipher=None

import os
base_path = os.getcwd()
icon_path = os.path.join(base_path, 'inti_logo4.ico')

a = Analysis(
    ['inti.py'],
    pathex=[],
    binaries=[],
    datas=[('inti_qt.ui','.'),('gong.ui','.'),
	('img_qt.ui','.'),('trame.ui','.'),('calc.ui','.'),
	('zoom.ui','.'),('grid.ui','.'),('profil_qt.ui','.'),
	('inti_logo.png','.'),('config_save.ui','.'), ('crop_box.ui','.'),('param.ui','.'),
	('shg_box.ui','.'),
	('matplotlib_cache/fontlist-v330.json', 'matplotlib_cache'),
	('matplotlib_cache/fontlist-v390.json', 'matplotlib_cache')],
    hiddenimports=[],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['PyQt5'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name='inti',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True, icon=icon_path,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)
coll = COLLECT(
    exe,
    a.binaries,
    a.datas,
    strip=False,
    upx=True,
    upx_exclude=[],
    name='inti',
)
