# -*- mode: python ; coding: utf-8 -*-
block_cipher=None

a = Analysis(
    ['inti_qt.py'],
    pathex=[],
    binaries=[],
    datas=[('inti_qt.ui','.'),('gong.ui','.'),
	('img_qt.ui','.'),('trame.ui','.'),('calc.ui','.'),
	('zoom.ui','.'),('grid.ui','.'),('profil_qt.ui','.'),
	('inti_logo.png','.'),('config_save.ui','.'), ('crop_box.ui','.'),('param.ui','.')],
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
    name='inti_qt',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    console=True, icon="C:\\Users\\valer\\codepy\\inti_qt\\inti_logo4.ico",
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
    name='inti_qt',
)
