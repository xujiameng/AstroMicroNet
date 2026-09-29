from pathlib import Path
import runpy
root=Path(__file__).resolve().parent
for name in ('plot_S1_S7.py','plot_updated_figures.py','plot_S9_quantitative.py'):
    runpy.run_path(str(root/name),run_name='__main__')
print('Exported S1-S8 and the quantitative-only S9 variant (panel A omitted).')
