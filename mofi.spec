a = Analysis(['run.py'],
             pathex=[],
             binaries=[],
             datas=[('mofi/config', 'config'),
                    ('mofi/data', 'data'),
                    ('mofi/docs', 'docs')],
             hiddenimports=['mofi.findmods',
                            'pandas._libs.tslibs.timedeltas'],
             hookspath=[],
             runtime_hooks=[],
             excludes=['jinja2'],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=None)

pyz = PYZ(a.pure,
          a.zipped_data,
          cipher=None)

exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='mofi',
          debug=False,
          strip=False,
          upx=True,
          console=True)

coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='mofi')
