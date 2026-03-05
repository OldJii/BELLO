"""
BELLO GUI  —  Cross-platform tkinter interface  (macOS / Windows / Linux)
GitHub-Dark theme  ·  embedded matplotlib  ·  Chinese/English i18n
"""

import platform
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import os
import time
from itertools import combinations_with_replacement

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import numpy as np

import BELLO_main
import Angle_Distribution_Function
import Radial_Pair_Distribution_Function
import Coordination_Heatmap
import vasp_converter

# ── Platform ──────────────────────────────────────────────────────────────────
IS_MAC = platform.system() == 'Darwin'
IS_WIN = platform.system() == 'Windows'


def _font_family():
    if IS_MAC: return 'Helvetica Neue'
    if IS_WIN: return 'Segoe UI'
    return 'DejaVu Sans'


_FF = _font_family()
F_XS  = (_FF,  8)
F_S   = (_FF,  9)
F_M   = (_FF, 10)
F_MB  = (_FF, 10, 'bold')
F_L   = (_FF, 11)
F_LB  = (_FF, 11, 'bold')
F_LOGO = (_FF, 21, 'bold')

# ── Color palette  (GitHub Dark) ──────────────────────────────────────────────
BG      = '#0d1117'
SURFACE = '#161b22'
PANEL   = '#1c2128'
RAISED  = '#21262d'
BORDER  = '#30363d'
BORDER2 = '#3d444d'
FG      = '#e6edf3'
FG2     = '#8b949e'
FG3     = '#656d76'
ACCENT  = '#2f81f7'
ACCENTL = '#58a6ff'
GREEN   = '#3fb950'
RED     = '#f85149'
ORANGE  = '#d29922'
INPUT   = '#0d1117'
SEL_BG  = '#1c3557'

# ── Matplotlib global style ───────────────────────────────────────────────────
plt.rcParams.update({
    'figure.facecolor':  'none',
    'axes.facecolor':    RAISED,
    'axes.edgecolor':    BORDER,
    'axes.labelcolor':   FG2,
    'axes.titlecolor':   FG,
    'text.color':        FG,
    'xtick.color':       FG2,
    'ytick.color':       FG2,
    'xtick.labelcolor':  FG2,
    'ytick.labelcolor':  FG2,
    'grid.color':        BORDER,
    'grid.linewidth':    0.6,
    'legend.facecolor':  RAISED,
    'legend.edgecolor':  BORDER,
    'legend.labelcolor': FG,
    'lines.linewidth':   1.8,
    'font.size':         9,
})


# ── Cancellation sentinel ─────────────────────────────────────────────────────
class _Cancelled(Exception):
    pass


# ── Thread-safe NavigationToolbar ─────────────────────────────────────────────
class _SafeToolbar(NavigationToolbar2Tk):
    def set_message(self, s):
        try:
            super().set_message(s)
        except RuntimeError:
            pass


# ── i18n ──────────────────────────────────────────────────────────────────────
_STRINGS = {
    'en': {
        'app_title':        'BELLO — Bond Element Lattice Locality Order',
        'logo':             'B.E.L.L.O',
        'subtitle':         'Local-order analysis for disordered systems',
        'lang_btn':         '中文',
        'sec_input':        'INPUT FILE',
        'browse':           'Browse',
        'sec_vasp':         'VASP → XYZ',
        'convert':          'Convert',
        'sec_params':       'PARAMETERS',
        'auto_thr':         'Automatic threshold',
        'threshold':        'Threshold',
        'tolerance':        'Tolerance',
        'frame_stride':     'Frame stride',
        'max_frames':       'Max frames  (0 = all)',
        'divide_angle':     'Split angle output  (memory-efficient)',
        'sec_cell':         'UNIT CELL  (Å)',
        'sec_calc':         'CALCULATION',
        'mode_bello':       'BELLO',
        'mode_rdf':         'RDF',
        'mode_angle':       'ADF',
        'mode_coord':       'Heatmap',
        'btn_calc':         'Run Calculation',
        'btn_stop':         'Stop',
        'btn_export':       'Export Excel',
        'ready':            'Ready',
        'running':          'Running…',
        'done':             'Done',
        'cancelled':        'Cancelled',
        'busy':             'A calculation is already in progress.',
        'error':            'Error',
        'plot_hint':        'Results will appear here after calculation.',
        'plot_tab':         'Plot',
        'dlg_rdf_title':    'RDF Parameters',
        'rdf_rmax':         'Max radius r',
        'rdf_dr':           'Step  Δr',
        'detected_elements':'Elements detected:  {}',
        'ok':               'OK',
        'cancel':           'Cancel',
        'err_vasp':         'Please select a valid VASP file.',
        'err_file':         'Please select a valid {} file.',
        'err_no_xyz':       'Please load an XYZ input file first.',
        'err_float':        '{} must be a number (got "{}")',
        'err_int':          '{} must be an integer (got "{}")',
        'err_min':          '{} must be ≥ {}',
        'conv_ok_title':    'Conversion Complete',
        'conv_ok_msg':      'Output: {}\nCell: {:.3f} × {:.3f} × {:.3f} Å\nFrames: {}',
        'export_ok':        'Exported: {}',
        'export_no_data':   'No data to export. Run a calculation first.',
        'sec_axis':         'AXIS RANGE',
        'x_min':            'X min',
        'x_max':            'X max',
        'btn_refresh':      'Refresh Plot',
        'sec_adf':          'ADF SETTINGS',
        'cutoff_auto':      'Auto cutoff (from RDF)',
        'cutoff_manual':    'Manual cutoff',
        'cutoff_value':     'Cutoff (Å)',
        'no_rdf_data':      'No RDF data available. Run RDF first to auto-detect cutoff, or use manual cutoff.',
        'heatmap_warn':     'Heatmap requires BELLO output data.\nPlease run BELLO calculation first.',
        'adf_no_cutoff':    'Cannot determine cutoff. Run RDF first or enter a manual cutoff value.',
    },
    'zh': {
        'app_title':        'BELLO — 键元素晶格局域有序度分析',
        'logo':             'B.E.L.L.O',
        'subtitle':         '无序/非晶体系局域有序度分析工具',
        'lang_btn':         'English',
        'sec_input':        '输入文件',
        'browse':           '浏览',
        'sec_vasp':         'VASP → XYZ',
        'convert':          '转换',
        'sec_params':       '参数设置',
        'auto_thr':         '自动确定阈值',
        'threshold':        '距离阈值',
        'tolerance':        '距离容差',
        'frame_stride':     '帧步长',
        'max_frames':       '最大帧数（0 = 全部）',
        'divide_angle':     '分段保存角度分布（节省内存）',
        'sec_cell':         '晶胞尺寸  (Å)',
        'sec_calc':         '计算',
        'mode_bello':       'BELLO',
        'mode_rdf':         'RDF',
        'mode_angle':       'ADF',
        'mode_coord':       '热力图',
        'btn_calc':         '开始计算',
        'btn_stop':         '停止计算',
        'btn_export':       '导出 Excel',
        'ready':            '就绪',
        'running':          '计算中…',
        'done':             '完成',
        'cancelled':        '已取消',
        'busy':             '计算任务正在进行中，请稍候。',
        'error':            '错误',
        'plot_hint':        '计算完成后，结果图表将显示在此处。',
        'plot_tab':         '图表',
        'dlg_rdf_title':    'RDF 参数',
        'rdf_rmax':         '最大半径 r',
        'rdf_dr':           '步长  Δr',
        'detected_elements':'检测到元素:  {}',
        'ok':               '确定',
        'cancel':           '取消',
        'err_vasp':         '请选择有效的 VASP 文件。',
        'err_file':         '请选择有效的 {} 文件。',
        'err_no_xyz':       '请先加载 XYZ 输入文件。',
        'err_float':        '{} 必须为数字（当前: "{}"）',
        'err_int':          '{} 必须为整数（当前: "{}"）',
        'err_min':          '{} 必须 ≥ {}',
        'conv_ok_title':    '转换成功',
        'conv_ok_msg':      '输出: {}\n晶胞: {:.3f} × {:.3f} × {:.3f} Å\n帧数: {}',
        'export_ok':        '已导出: {}',
        'export_no_data':   '没有可导出的数据，请先运行计算。',
        'sec_axis':         '坐标轴范围',
        'x_min':            'X 最小值',
        'x_max':            'X 最大值',
        'btn_refresh':      '刷新图表',
        'sec_adf':          'ADF 设置',
        'cutoff_auto':      '自动截距（来自 RDF）',
        'cutoff_manual':    '手动截距',
        'cutoff_value':     '截距 (Å)',
        'no_rdf_data':      '无 RDF 数据。请先运行 RDF 以自动检测截距，或使用手动截距。',
        'heatmap_warn':     '热力图需要 BELLO 输出数据。\n请先运行 BELLO 计算。',
        'adf_no_cutoff':    '无法确定截距。请先运行 RDF 或手动输入截距值。',
    },
}


# ═══════════════════════════════════════════════════════════════════════════════
# Main application
# ═══════════════════════════════════════════════════════════════════════════════
class BelloApp(tk.Tk):

    def __init__(self):
        super().__init__()

        self._lang         = 'en'
        self._t            = _STRINGS['en']
        self._busy         = False
        self._cancel_event = threading.Event()
        self._prog_ts      = 0.0
        self._prog_stage   = ''
        self._prog_pct     = -1

        # Cached computation results for cross-module linkage
        self._rdf_cache = {}       # {(e1,e2): {'r': array, 'rdf': array, 'first_min': float|None}}
        self._last_figures = []    # last displayed figures
        self._last_calc_type = ''  # 'BELLO', 'RDF', 'ANGLE', 'COORD'
        self._last_calc_label = '' # e.g. 'Ge-Se' for export naming
        self._last_export_data = {}  # {'type': ..., 'data': {name: (x_arr, y_arr)}}

        self.title(self._t['app_title'])
        self.configure(bg=BG)
        self.minsize(1020, 680)
        self.geometry('1280x820')

        if IS_MAC:
            self.createcommand('tk::mac::Quit', self.destroy)

        self._apply_theme()
        self._build_ui()

    # ── i18n ──────────────────────────────────────────────────────────────────
    def t(self, key):
        return self._t.get(key, key)

    def _switch_lang(self):
        self._lang = 'zh' if self._lang == 'en' else 'en'
        self._t = _STRINGS[self._lang]
        self.title(self.t('app_title'))
        self._rebuild_all()

    # ── ttk theme ─────────────────────────────────────────────────────────────
    def _apply_theme(self):
        s = ttk.Style(self)
        s.theme_use('clam')

        s.configure('.',              background=BG,     foreground=FG,  borderwidth=0, font=F_M)
        s.configure('TFrame',         background=BG)
        s.configure('Panel.TFrame',   background=PANEL)
        s.configure('Raised.TFrame',  background=RAISED)
        s.configure('Surface.TFrame', background=SURFACE)

        s.configure('TLabel',        background=BG,     foreground=FG,  font=F_M)
        s.configure('Panel.TLabel',  background=PANEL,  foreground=FG,  font=F_M)
        s.configure('Muted.TLabel',  background=PANEL,  foreground=FG2, font=F_S)
        s.configure('Hint.TLabel',   background=RAISED, foreground=FG2, font=F_L)

        s.configure('TEntry',
                    fieldbackground=INPUT, foreground=FG,
                    insertcolor=FG, relief='flat', borderwidth=0,
                    padding=(6, 4), font=F_M)
        s.map('TEntry',
              fieldbackground=[('focus', '#0b1d3a')],
              lightcolor=[('focus', ACCENT)],
              darkcolor=[('focus', ACCENT)])

        for w in ('TCheckbutton', 'TRadiobutton'):
            s.configure(w, background=PANEL, foreground=FG, font=F_M)
            s.map(w,
                  background=[('active', PANEL), ('!active', PANEL)],
                  foreground=[('disabled', FG3)])

        s.configure('Primary.TButton',
                    background=ACCENT, foreground='#ffffff',
                    font=F_MB, padding=(0, 10), relief='flat', borderwidth=0)
        s.map('Primary.TButton',
              background=[('active', ACCENTL), ('disabled', RAISED)],
              foreground=[('disabled', FG3)])

        s.configure('Danger.TButton',
                    background=RED, foreground='#ffffff',
                    font=F_MB, padding=(0, 10), relief='flat', borderwidth=0)
        s.map('Danger.TButton',
              background=[('active', '#ff6b6b'), ('disabled', RAISED)],
              foreground=[('disabled', FG3)])

        s.configure('Secondary.TButton',
                    background=RAISED, foreground=FG,
                    font=F_S, padding=(8, 5), relief='flat', borderwidth=1)
        s.map('Secondary.TButton',
              background=[('active', BORDER), ('disabled', PANEL)],
              foreground=[('disabled', FG3)])

        s.configure('Export.TButton',
                    background='#238636', foreground='#ffffff',
                    font=F_S, padding=(8, 5), relief='flat', borderwidth=0)
        s.map('Export.TButton',
              background=[('active', '#2ea043'), ('disabled', PANEL)],
              foreground=[('disabled', FG3)])

        s.configure('Lang.TButton',
                    background=SURFACE, foreground=ACCENTL,
                    font=F_XS, padding=(9, 4), relief='flat')
        s.map('Lang.TButton', background=[('active', PANEL)])

        s.configure('TNotebook',
                    background=BG, borderwidth=0, tabmargins=(0, 0, 0, 0))
        s.configure('TNotebook.Tab',
                    background=SURFACE, foreground=FG2,
                    font=F_S, padding=(18, 8), borderwidth=0)
        s.map('TNotebook.Tab',
              background=[('selected', RAISED)],
              foreground=[('selected', FG)])

        s.configure('Thin.Horizontal.TProgressbar',
                    troughcolor=RAISED, background=ACCENT,
                    borderwidth=0, thickness=2)

        s.configure('TScrollbar',
                    background=RAISED, troughcolor=PANEL,
                    borderwidth=0, arrowsize=10)
        s.map('TScrollbar', background=[('active', BORDER)])

        s.configure('TCombobox',
                    fieldbackground=INPUT, foreground=FG,
                    background=RAISED, selectbackground=SEL_BG,
                    selectforeground=FG, font=F_M)
        s.map('TCombobox',
              fieldbackground=[('readonly', INPUT)],
              selectbackground=[('readonly', SEL_BG)])

    # ── build / rebuild ────────────────────────────────────────────────────────
    def _rebuild_all(self):
        for w in self.winfo_children():
            w.destroy()
        self._apply_theme()
        self._build_ui()

    def _build_ui(self):
        self._build_topbar()
        self._build_main()
        self._build_statusbar()

    # ── top bar ────────────────────────────────────────────────────────────────
    def _build_topbar(self):
        bar = tk.Frame(self, bg=SURFACE, height=54)
        bar.pack(fill='x')
        bar.pack_propagate(False)

        inner = tk.Frame(bar, bg=SURFACE)
        inner.pack(fill='both', expand=True, padx=20)

        tk.Label(inner, text=self.t('logo'),
                 bg=SURFACE, fg=FG, font=F_LOGO).pack(side='left', pady=10)
        tk.Label(inner, text='·',
                 bg=SURFACE, fg=FG3, font=(_FF, 16)).pack(side='left', padx=8, pady=10)
        tk.Label(inner, text=self.t('subtitle'),
                 bg=SURFACE, fg=FG2, font=F_S).pack(side='left', pady=(18, 0))

        self.btn_lang = ttk.Button(inner, text=self.t('lang_btn'),
                                   style='Lang.TButton',
                                   command=self._switch_lang)
        self.btn_lang.pack(side='right', pady=13)

        tk.Frame(self, bg=BORDER, height=1).pack(fill='x')

    # ── main layout ───────────────────────────────────────────────────────────
    def _build_main(self):
        main = tk.Frame(self, bg=BG)
        main.pack(fill='both', expand=True)

        self._build_panel(main)
        tk.Frame(main, bg=BORDER, width=1).pack(side='left', fill='y')
        self._build_plot_area(main)

    # ── left control panel ────────────────────────────────────────────────────
    def _build_panel(self, parent):
        outer = tk.Frame(parent, bg=PANEL, width=330)
        outer.pack(side='left', fill='y')
        outer.pack_propagate(False)

        canvas = tk.Canvas(outer, bg=PANEL, highlightthickness=0, width=328)
        sb = ttk.Scrollbar(outer, orient='vertical', command=canvas.yview)
        self._scroll_inner = tk.Frame(canvas, bg=PANEL)

        self._scroll_inner.bind(
            '<Configure>',
            lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=self._scroll_inner, anchor='nw')
        canvas.configure(yscrollcommand=sb.set)

        canvas.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')

        self._panel_canvas = canvas
        self._bind_scroll(canvas)

        self._build_input_section(self._scroll_inner)
        self._build_vasp_section(self._scroll_inner)
        self._build_params_section(self._scroll_inner)
        self._build_cell_section(self._scroll_inner)
        self._build_adf_section(self._scroll_inner)
        self._build_axis_section(self._scroll_inner)
        self._build_calc_section(self._scroll_inner)

    def _bind_scroll(self, canvas):
        if IS_MAC:
            canvas.bind_all('<MouseWheel>',
                            lambda e: canvas.yview_scroll(-e.delta, 'units'))
        elif IS_WIN:
            canvas.bind_all('<MouseWheel>',
                            lambda e: canvas.yview_scroll(int(-e.delta / 120), 'units'))
        else:
            canvas.bind_all('<Button-4>', lambda e: canvas.yview_scroll(-3, 'units'))
            canvas.bind_all('<Button-5>', lambda e: canvas.yview_scroll(3,  'units'))

    # ── section factory ───────────────────────────────────────────────────────
    def _section(self, parent, title_key):
        sec = tk.Frame(parent, bg=PANEL)
        sec.pack(fill='x')

        hdr = tk.Frame(sec, bg=PANEL)
        hdr.pack(fill='x', padx=16, pady=(20, 6))
        tk.Label(hdr, text=self.t(title_key),
                 bg=PANEL, fg=ACCENTL, font=F_XS).pack(side='left')
        tk.Frame(hdr, bg=BORDER, height=1).pack(side='right', fill='x', expand=True,
                                                 padx=(10, 0), pady=4)
        body = tk.Frame(sec, bg=PANEL)
        body.pack(fill='x', padx=16, pady=(0, 14))
        return body

    def _input_row(self, parent, label, default='', width=9):
        row = tk.Frame(parent, bg=PANEL)
        row.pack(fill='x', pady=3)
        tk.Label(row, text=label, bg=PANEL, fg=FG2, font=F_S).pack(side='left')
        var = tk.StringVar(value=default)
        ef = tk.Frame(row, bg=BORDER2, padx=1, pady=1)
        ef.pack(side='right')
        ttk.Entry(ef, textvariable=var, width=width, font=F_M).pack()
        return var

    # ── INPUT FILE section ────────────────────────────────────────────────────
    def _build_input_section(self, parent):
        body = self._section(parent, 'sec_input')
        row = tk.Frame(body, bg=PANEL)
        row.pack(fill='x')
        self.url_var = tk.StringVar()
        ef = tk.Frame(row, bg=BORDER2, padx=1, pady=1)
        ef.pack(side='left', fill='x', expand=True)
        ttk.Entry(ef, textvariable=self.url_var, font=F_S).pack(fill='x')
        self.btn_browse_xyz = ttk.Button(row, text=self.t('browse'),
                                         style='Secondary.TButton',
                                         command=self._browse_xyz)
        self.btn_browse_xyz.pack(side='right', padx=(6, 0))

    def _browse_xyz(self):
        p = filedialog.askopenfilename(title='XYZ')
        if p:
            self.url_var.set(p)

    # ── VASP → XYZ section ────────────────────────────────────────────────────
    def _build_vasp_section(self, parent):
        body = self._section(parent, 'sec_vasp')
        row = tk.Frame(body, bg=PANEL)
        row.pack(fill='x')
        self.vasp_var = tk.StringVar()
        ef = tk.Frame(row, bg=BORDER2, padx=1, pady=1)
        ef.pack(side='left', fill='x', expand=True)
        ttk.Entry(ef, textvariable=self.vasp_var, font=F_S).pack(fill='x')
        self.btn_browse_vasp = ttk.Button(row, text=self.t('browse'),
                                          style='Secondary.TButton',
                                          command=self._browse_vasp)
        self.btn_browse_vasp.pack(side='right', padx=(6, 0))
        self.btn_convert_vasp = ttk.Button(body, text=self.t('convert'),
                                           style='Secondary.TButton',
                                           command=self._convert_vasp)
        self.btn_convert_vasp.pack(fill='x', pady=(8, 0))

    def _browse_vasp(self):
        p = filedialog.askopenfilename(title='VASP')
        if p:
            self.vasp_var.set(p)

    def _convert_vasp(self):
        vp = self.vasp_var.get().strip()
        if not vp or not os.path.isfile(vp):
            messagebox.showerror(self.t('error'), self.t('err_vasp'))
            return
        try:
            xyz, cx, cy, cz, nframes = vasp_converter.convert_vasp(vp)
            self.url_var.set(xyz)
            self.cell_x.set(f'{cx:.5f}')
            self.cell_y.set(f'{cy:.5f}')
            self.cell_z.set(f'{cz:.5f}')
            messagebox.showinfo(
                self.t('conv_ok_title'),
                self.t('conv_ok_msg').format(xyz, cx, cy, cz, nframes))
        except Exception as e:
            messagebox.showerror(self.t('error'), str(e))

    # ── PARAMETERS section ────────────────────────────────────────────────────
    def _build_params_section(self, parent):
        body = self._section(parent, 'sec_params')
        self.auto_var = tk.BooleanVar(value=False)
        self.chk_auto = ttk.Checkbutton(body, text=self.t('auto_thr'),
                                        variable=self.auto_var,
                                        command=self._toggle_auto)
        self.chk_auto.pack(anchor='w', pady=(0, 6))
        self.trh_var          = self._input_row(body, self.t('threshold'),    '3.0')
        self.tlr_var          = self._input_row(body, self.t('tolerance'),    '0.5')
        self.frame_stride_var = self._input_row(body, self.t('frame_stride'), '1')
        self.max_frames_var   = self._input_row(body, self.t('max_frames'),   '0')
        self.sep_ang_var = tk.BooleanVar(value=False)
        self.chk_sep_angle = ttk.Checkbutton(body, text=self.t('divide_angle'),
                                             variable=self.sep_ang_var)
        self.chk_sep_angle.pack(anchor='w', pady=(8, 0))

    def _toggle_auto(self):
        st = 'disabled' if self.auto_var.get() else 'normal'
        for w in self._all_widgets():
            if isinstance(w, ttk.Entry):
                try:
                    if str(w.cget('textvariable')) in (str(self.trh_var), str(self.tlr_var)):
                        w.configure(state=st)
                except Exception:
                    pass

    def _all_widgets(self):
        stack = list(self.winfo_children())
        while stack:
            w = stack.pop()
            yield w
            stack.extend(w.winfo_children())

    # ── UNIT CELL section ─────────────────────────────────────────────────────
    def _build_cell_section(self, parent):
        body = self._section(parent, 'sec_cell')
        self.cell_x = self._input_row(body, 'X', '20.00000')
        self.cell_y = self._input_row(body, 'Y', '20.00000')
        self.cell_z = self._input_row(body, 'Z', '20.00000')

    # ── ADF SETTINGS section ──────────────────────────────────────────────────
    def _build_adf_section(self, parent):
        body = self._section(parent, 'sec_adf')

        self.cutoff_mode = tk.StringVar(value='auto')
        ttk.Radiobutton(body, text=self.t('cutoff_auto'),
                        variable=self.cutoff_mode, value='auto',
                        command=self._toggle_cutoff_mode).pack(anchor='w')
        ttk.Radiobutton(body, text=self.t('cutoff_manual'),
                        variable=self.cutoff_mode, value='manual',
                        command=self._toggle_cutoff_mode).pack(anchor='w', pady=(2, 6))
        self.cutoff_var = self._input_row(body, self.t('cutoff_value'), '', width=9)

        self._cutoff_info = tk.Label(body, text='', bg=PANEL, fg=FG3, font=F_XS,
                                     wraplength=260, justify='left')
        self._cutoff_info.pack(anchor='w', pady=(4, 0))
        self._update_cutoff_info()

    def _toggle_cutoff_mode(self):
        is_manual = (self.cutoff_mode.get() == 'manual')
        for w in self._all_widgets():
            if isinstance(w, ttk.Entry):
                try:
                    if str(w.cget('textvariable')) == str(self.cutoff_var):
                        w.configure(state='normal' if is_manual else 'disabled')
                except Exception:
                    pass
        self._update_cutoff_info()

    def _update_cutoff_info(self):
        if not hasattr(self, '_cutoff_info'):
            return
        if self.cutoff_mode.get() == 'auto':
            if self._rdf_cache:
                lines = []
                for (e1, e2), data in self._rdf_cache.items():
                    fm = data.get('first_min')
                    lines.append(f'{e1}–{e2}: {fm:.2f} Å' if fm else f'{e1}–{e2}: N/A')
                self._cutoff_info.config(text='RDF cutoffs:\n' + '\n'.join(lines), fg=GREEN)
            else:
                self._cutoff_info.config(text=self.t('no_rdf_data'), fg=ORANGE)
        else:
            self._cutoff_info.config(text='', fg=FG3)

    # ── AXIS RANGE section ────────────────────────────────────────────────────
    def _build_axis_section(self, parent):
        body = self._section(parent, 'sec_axis')
        self.xmin_var = self._input_row(body, self.t('x_min'), '', width=7)
        self.xmax_var = self._input_row(body, self.t('x_max'), '', width=7)
        btn_row = tk.Frame(body, bg=PANEL)
        btn_row.pack(fill='x', pady=(6, 0))
        self.btn_refresh = ttk.Button(btn_row, text=self.t('btn_refresh'),
                                      style='Secondary.TButton',
                                      command=self._refresh_axis)
        self.btn_refresh.pack(fill='x')

    def _refresh_axis(self):
        if not self._last_figures:
            return
        xmin_s = self.xmin_var.get().strip()
        xmax_s = self.xmax_var.get().strip()
        xmin = float(xmin_s) if xmin_s else None
        xmax = float(xmax_s) if xmax_s else None
        for fig in self._last_figures:
            for ax in fig.get_axes():
                if xmin is not None or xmax is not None:
                    cur = ax.get_xlim()
                    ax.set_xlim(xmin if xmin is not None else cur[0],
                                xmax if xmax is not None else cur[1])
        self._show_figures(self._last_figures, keep_cache=True)

    # ── CALCULATION section ────────────────────────────────────────────────────
    def _build_calc_section(self, parent):
        body = self._section(parent, 'sec_calc')

        self.calc_mode = tk.StringVar(value='BELLO')
        mode_defs = [
            ('mode_bello', 'BELLO'),
            ('mode_rdf',   'RDF'),
            ('mode_angle', 'ANGLE'),
            ('mode_coord', 'COORD'),
        ]

        seg_outer = tk.Frame(body, bg=RAISED)
        seg_outer.pack(fill='x', pady=(0, 14))
        tk.Frame(seg_outer, bg=BORDER, height=1).pack(fill='x')

        seg = tk.Frame(seg_outer, bg=RAISED)
        seg.pack(fill='x')

        self._mode_btns = {}
        for col, (key, val) in enumerate(mode_defs):
            is_sel = (val == 'BELLO')
            lbl = tk.Label(seg,
                           text=self.t(key),
                           bg=SEL_BG if is_sel else RAISED,
                           fg=ACCENTL if is_sel else FG2,
                           font=F_S, cursor='hand2',
                           padx=0, pady=8, anchor='center')
            lbl.grid(row=0, column=col, sticky='ew')
            seg.columnconfigure(col, weight=1)
            lbl.bind('<Button-1>', lambda e, v=val: self._select_mode(v))
            lbl.bind('<Enter>',    lambda e, w=lbl, v=val:
                     w.config(bg=BORDER if v != self.calc_mode.get() else SEL_BG))
            lbl.bind('<Leave>',    lambda e, w=lbl, v=val:
                     w.config(bg=SEL_BG if v == self.calc_mode.get() else RAISED))
            self._mode_btns[val] = lbl

        tk.Frame(seg_outer, bg=BORDER, height=1).pack(fill='x')

        # Run + Stop + Export buttons
        self._btn_slot = tk.Frame(body, bg=PANEL)
        self._btn_slot.pack(fill='x')

        self.btn_calc = ttk.Button(self._btn_slot, text=self.t('btn_calc'),
                                   style='Primary.TButton',
                                   command=self._on_calculate)
        self.btn_calc.pack(fill='x')

        self.btn_stop = ttk.Button(self._btn_slot, text=self.t('btn_stop'),
                                   style='Danger.TButton',
                                   command=self._stop_calculation)

        self.btn_export = ttk.Button(body, text=self.t('btn_export'),
                                     style='Export.TButton',
                                     command=self._export_excel)
        self.btn_export.pack(fill='x', pady=(8, 0))

    def _select_mode(self, val):
        if self._busy:
            return
        self.calc_mode.set(val)
        for v, lbl in self._mode_btns.items():
            lbl.config(bg=SEL_BG if v == val else RAISED,
                       fg=ACCENTL if v == val else FG2)

    # ── Plot area ─────────────────────────────────────────────────────────────
    def _build_plot_area(self, parent):
        self.plot_notebook = ttk.Notebook(parent)
        self.plot_notebook.pack(side='left', fill='both', expand=True)

        ph = tk.Frame(self.plot_notebook, bg=RAISED)
        self.plot_notebook.add(ph, text=f'  {self.t("plot_tab")}  ')
        tk.Label(ph, text=self.t('plot_hint'),
                 bg=RAISED, fg=FG2, font=F_L).place(relx=0.5, rely=0.5, anchor='center')

    def _show_figures(self, figs, keep_cache=False):
        for tab_id in self.plot_notebook.tabs():
            self.plot_notebook.forget(tab_id)
        for idx, fig in enumerate(figs):
            frame = tk.Frame(self.plot_notebook, bg=RAISED)
            title = f'  {self.t("plot_tab")} {idx + 1}  '
            for ax in fig.get_axes():
                t = ax.get_title()
                if t:
                    title = f'  {t[:30]}  '
                    break
            self.plot_notebook.add(frame, text=title)
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.draw()
            tb_frame = tk.Frame(frame, bg=RAISED)
            tb_frame.pack(side='bottom', fill='x')
            toolbar = _SafeToolbar(canvas, tb_frame)
            try:
                toolbar.config(bg=RAISED)
                for child in toolbar.winfo_children():
                    child.config(bg=RAISED)
            except Exception:
                pass
            canvas.get_tk_widget().config(bg=RAISED, highlightthickness=0)
            canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
        if figs:
            self.plot_notebook.select(0)
        if not keep_cache:
            self._last_figures = list(figs)

    # ── Status bar ────────────────────────────────────────────────────────────
    def _build_statusbar(self):
        tk.Frame(self, bg=BORDER, height=1).pack(fill='x')
        bar = tk.Frame(self, bg=SURFACE, height=30)
        bar.pack(fill='x')
        bar.pack_propagate(False)

        inner = tk.Frame(bar, bg=SURFACE)
        inner.pack(fill='both', expand=True, padx=16)

        self._dot = tk.Label(inner, text='●', bg=SURFACE, fg=GREEN, font=F_XS)
        self._dot.pack(side='left', pady=8)

        self.status_var = tk.StringVar(value=self.t('ready'))
        tk.Label(inner, textvariable=self.status_var,
                 bg=SURFACE, fg=FG2, font=F_S).pack(side='left', padx=(5, 0), pady=8)

        self.progress = ttk.Progressbar(inner, mode='indeterminate',
                                        style='Thin.Horizontal.TProgressbar',
                                        length=220)
        self.progress.pack(side='right', pady=14)

    # ── Validation ────────────────────────────────────────────────────────────
    def _vfloat(self, value, name='value'):
        try:
            return float(value)
        except (ValueError, TypeError):
            raise ValueError(self.t('err_float').format(name, value))

    def _vint(self, value, name='value', min_value=None):
        try:
            ivalue = int(value)
        except (ValueError, TypeError):
            raise ValueError(self.t('err_int').format(name, value))
        if min_value is not None and ivalue < min_value:
            raise ValueError(self.t('err_min').format(name, min_value))
        return ivalue

    def _vfile(self, path, label='file'):
        if not path or not os.path.isfile(path):
            raise ValueError(self.t('err_file').format(label))

    # ── Element auto-detection ────────────────────────────────────────────────
    @staticmethod
    def _detect_elements(xyz_path):
        elements = set()
        with open(xyz_path, 'r') as f:
            lines = f.readlines()
        i = 0
        while i < len(lines):
            line = lines[i].strip()
            if not line:
                i += 1
                continue
            try:
                n = int(line)
                i += 2
                for j in range(min(n, len(lines) - i)):
                    parts = lines[i + j].split()
                    if parts:
                        elements.add(parts[0])
                break
            except ValueError:
                i += 1
        return sorted(elements)

    # ── Busy state ────────────────────────────────────────────────────────────
    def _set_busy(self, busy):
        self._busy = busy
        state = 'disabled' if busy else 'normal'

        for w in [
            getattr(self, 'btn_lang',         None),
            getattr(self, 'btn_browse_xyz',   None),
            getattr(self, 'btn_browse_vasp',  None),
            getattr(self, 'btn_convert_vasp', None),
            getattr(self, 'chk_auto',         None),
            getattr(self, 'chk_sep_angle',    None),
            getattr(self, 'btn_export',       None),
            getattr(self, 'btn_refresh',      None),
        ]:
            if w is not None and w.winfo_exists():
                try:
                    w.configure(state=state)
                except tk.TclError:
                    pass

        if hasattr(self, 'btn_calc') and self.btn_calc.winfo_exists():
            if busy:
                self.btn_calc.pack_forget()
                self.btn_stop.pack(fill='x')
            else:
                self.btn_stop.pack_forget()
                self.btn_calc.pack(fill='x')

        cur = self.calc_mode.get()
        for val, lbl in getattr(self, '_mode_btns', {}).items():
            if lbl.winfo_exists():
                if busy:
                    lbl.config(fg=FG3, cursor='')
                else:
                    lbl.config(fg=ACCENTL if val == cur else FG2, cursor='hand2')

        if hasattr(self, '_dot') and self._dot.winfo_exists():
            self._dot.config(fg=RED if busy else GREEN)

        if not busy:
            self._toggle_auto()
            self._toggle_cutoff_mode()

    # ── Stop calculation ──────────────────────────────────────────────────────
    def _stop_calculation(self):
        self._cancel_event.set()
        if hasattr(self, 'btn_stop') and self.btn_stop.winfo_exists():
            self.btn_stop.configure(state='disabled')
        self.status_var.set(self.t('cancelled') + '…')

    # ── Calculation dispatch ──────────────────────────────────────────────────
    def _on_calculate(self):
        if self._busy:
            self.status_var.set(self.t('busy'))
            return
        mode = self.calc_mode.get()
        try:
            {'BELLO': self._run_bello,
             'RDF':   self._run_rdf,
             'ANGLE': self._run_angle,
             'COORD': self._run_coord}[mode]()
        except ValueError as e:
            messagebox.showerror(self.t('error'), str(e))
        except Exception as e:
            messagebox.showerror(self.t('error'), str(e))

    def _threaded(self, fn, args=()):
        self._set_busy(True)
        self._cancel_event.clear()
        self.progress.start(12)
        self.status_var.set(self.t('running'))
        self._prog_ts    = 0.0
        self._prog_stage = ''
        self._prog_pct   = -1

        def _worker():
            try:
                result = fn(*args)
                self.after(0, self._done, result)
            except _Cancelled:
                self.after(0, self._task_cancelled)
            except Exception as exc:
                self.after(0, self._fail, exc)

        threading.Thread(target=_worker, daemon=True).start()

    def _done(self, figures):
        self.progress.stop()
        self._cancel_event.clear()
        self._set_busy(False)
        self.status_var.set(self.t('done'))
        if figures:
            self._show_figures(figures)

    def _task_cancelled(self):
        self.progress.stop()
        self._cancel_event.clear()
        self._set_busy(False)
        self.status_var.set(self.t('cancelled'))

    def _fail(self, exc):
        self.progress.stop()
        self._cancel_event.clear()
        self._set_busy(False)
        self.status_var.set(self.t('error'))
        messagebox.showerror(self.t('error'), str(exc))

    def _pcb(self, stage, cur, total):
        if self._cancel_event.is_set():
            raise _Cancelled()
        pct = int(cur / max(total, 1) * 100)
        now = time.monotonic()
        if stage == self._prog_stage and pct < 100:
            if pct == self._prog_pct:
                return
            if now - self._prog_ts < 0.1:
                return
        self._prog_stage = stage
        self._prog_pct   = pct
        self._prog_ts    = now
        self.after(0, lambda: self.status_var.set(f'{stage}  {pct}%'))

    # ── BELLO ─────────────────────────────────────────────────────────────────
    def _run_bello(self):
        url = self.url_var.get().strip()
        self._vfile(url, 'XYZ')
        cx = self._vfloat(self.cell_x.get(), 'X')
        cy = self._vfloat(self.cell_y.get(), 'Y')
        cz = self._vfloat(self.cell_z.get(), 'Z')
        fs = self._vint(self.frame_stride_var.get(), self.t('frame_stride'), min_value=1)
        mf = self._vint(self.max_frames_var.get(),   self.t('max_frames'),   min_value=0)
        auto = self.auto_var.get()
        trh  = self.trh_var.get() if not auto else '2'
        tlr  = self.tlr_var.get() if not auto else '1'
        self._last_calc_type = 'BELLO'
        self._last_calc_label = 'BELLO'
        self._threaded(BELLO_main.BELLO,
                       (url, cx, cy, cz, auto, trh, tlr,
                        self.sep_ang_var.get(), self._pcb, fs, mf))

    # ── RDF  (all element pairs, auto-detected) ───────────────────────────────
    def _run_rdf(self):
        url = self.url_var.get().strip()
        self._vfile(url, 'XYZ')
        cx = self._vfloat(self.cell_x.get(), 'X')
        cy = self._vfloat(self.cell_y.get(), 'Y')
        cz = self._vfloat(self.cell_z.get(), 'Z')
        fs = self._vint(self.frame_stride_var.get(), self.t('frame_stride'), min_value=1)
        mf = self._vint(self.max_frames_var.get(),   self.t('max_frames'),   min_value=0)

        try:
            elements = self._detect_elements(url)
        except Exception:
            elements = []
        if not elements:
            raise ValueError(self.t('err_no_xyz'))

        dlg = _RdfParamsDialog(self, self._t, elements)
        self.wait_window(dlg)
        if dlg.result is None:
            return
        dr, rmax = dlg.result
        self._last_calc_type = 'RDF'
        self._last_calc_label = 'all'
        self._threaded(self._do_rdf_pairs,
                       (url, cx, cy, cz, dr, rmax, elements, fs, mf))

    def _do_rdf_pairs(self, url, cx, cy, cz, dr, rmax, elements, fs, mf):
        pairs = list(combinations_with_replacement(elements, 2))
        all_figs = []
        export_data = {}
        for idx, (e1, e2) in enumerate(pairs):
            def pair_pcb(stage, cur, total, _e1=e1, _e2=e2):
                self._pcb(f'{_e1}–{_e2}  {stage}', cur, total)
            try:
                figs, r_arr, rdf_arr, first_min = Radial_Pair_Distribution_Function.RDF(
                    url, cx, cy, cz, dr, rmax, e1, e2, pair_pcb, fs, mf)
                self._rdf_cache[(e1, e2)] = {
                    'r': r_arr, 'rdf': rdf_arr, 'first_min': first_min
                }
                if e1 != e2:
                    self._rdf_cache[(e2, e1)] = self._rdf_cache[(e1, e2)]
                for fig in figs:
                    for ax in fig.get_axes():
                        ax.set_title(f'g(r)  {e1}–{e2}', fontsize=11, color=FG)
                all_figs.extend(figs)
                export_data[f'{e1}-{e2}'] = (r_arr, rdf_arr)
            except _Cancelled:
                raise
            except Exception as exc:
                print(f"RDF pair {e1}-{e2} skipped: {exc}")

        self._last_export_data = {'type': 'RDF', 'data': export_data}
        self.after(0, self._update_cutoff_info)

        if not all_figs:
            raise ValueError("No RDF data computed for any element pair.")
        return all_figs

    # ── ADF  (direct computation from XYZ using cutoff) ──────────────────────
    def _run_angle(self):
        url = self.url_var.get().strip()
        self._vfile(url, 'XYZ')
        cx = self._vfloat(self.cell_x.get(), 'X')
        cy = self._vfloat(self.cell_y.get(), 'Y')
        cz = self._vfloat(self.cell_z.get(), 'Z')
        fs = self._vint(self.frame_stride_var.get(), self.t('frame_stride'), min_value=1)
        mf = self._vint(self.max_frames_var.get(),   self.t('max_frames'),   min_value=0)

        try:
            elements = self._detect_elements(url)
        except Exception:
            elements = []
        if not elements:
            raise ValueError(self.t('err_no_xyz'))

        is_auto = (self.cutoff_mode.get() == 'auto')

        if is_auto and not self._rdf_cache:
            raise ValueError(self.t('adf_no_cutoff'))

        self._last_calc_type = 'ANGLE'
        self._last_calc_label = 'all'
        self._threaded(self._do_adf_pairs,
                       (url, cx, cy, cz, elements, is_auto, fs, mf))

    def _do_adf_pairs(self, url, cx, cy, cz, elements, is_auto, fs, mf):
        all_figs = []
        export_data = {}

        triplets = []
        for center in elements:
            ligand_pairs = list(combinations_with_replacement(elements, 2))
            for (l1, l2) in ligand_pairs:
                triplets.append((l1, center, l2))

        for (l1, center, l2) in triplets:
            if is_auto:
                cutoff = self._get_auto_cutoff(center, l1, l2)
                if cutoff is None:
                    print(f"ADF {l1}-{center}-{l2}: no cutoff available, skipping")
                    continue
            else:
                cutoff_s = self.cutoff_var.get().strip()
                if not cutoff_s:
                    raise ValueError(self.t('adf_no_cutoff'))
                cutoff = float(cutoff_s)

            triplet_name = f'{l1}-{center}-{l2}'
            def trip_pcb(stage, cur, total, _name=triplet_name):
                self._pcb(f'{_name}  {stage}', cur, total)

            try:
                angles, fig = Angle_Distribution_Function.compute_adf(
                    url, cx, cy, cz, center, l1, l2,
                    cutoff, trip_pcb, fs, mf)
                if angles.size > 0:
                    for ax in fig.get_axes():
                        ax.set_title(f'ADF {triplet_name} (cutoff={cutoff:.2f}Å)',
                                     fontsize=10, color=FG)
                    all_figs.append(fig)
                    export_data[triplet_name] = (angles, None)
                else:
                    plt.close(fig)
            except _Cancelled:
                raise
            except Exception as exc:
                print(f"ADF {triplet_name} skipped: {exc}")

        self._last_export_data = {'type': 'ADF', 'data': export_data}
        if not all_figs:
            raise ValueError("No ADF data computed for any triplet.")
        return all_figs

    def _get_auto_cutoff(self, center, l1, l2):
        """Find the appropriate cutoff from RDF cache for a triplet l1-center-l2.

        For the triplet Se-Ge-Se, we need the Ge-Se RDF first minimum.
        We look up (center, l1) and (center, l2) and take the max.
        """
        cutoffs = []
        for ligand in set([l1, l2]):
            key1 = (center, ligand)
            key2 = (ligand, center)
            data = self._rdf_cache.get(key1) or self._rdf_cache.get(key2)
            if data and data.get('first_min') is not None:
                cutoffs.append(data['first_min'])
        return max(cutoffs) if cutoffs else None

    # ── Coordination Heatmap  (auto-detected elements) ────────────────────────
    def _run_coord(self):
        url = self.url_var.get().strip()
        if not url or not os.path.isfile(url):
            raise ValueError(self.t('err_no_xyz'))

        if not os.path.isfile('output-human-readable-coords.txt'):
            messagebox.showwarning(self.t('error'), self.t('heatmap_warn'))
            return

        try:
            elements = self._detect_elements(url)
        except Exception:
            elements = []
        if not elements:
            raise ValueError(self.t('err_no_xyz'))
        self._last_calc_type = 'COORD'
        self._last_calc_label = 'all'
        self._threaded(Coordination_Heatmap.coordination_heatmap, (elements,))

    # ── Excel Export ──────────────────────────────────────────────────────────
    def _export_excel(self):
        if not self._last_export_data or not self._last_export_data.get('data'):
            messagebox.showinfo(self.t('error'), self.t('export_no_data'))
            return

        try:
            import openpyxl
        except ImportError:
            messagebox.showerror(self.t('error'),
                                 'openpyxl is required for Excel export.\n'
                                 'Install: pip install openpyxl')
            return

        calc_type = self._last_export_data['type']
        data_dict = self._last_export_data['data']
        exported = []

        for name, (arr_x, arr_y) in data_dict.items():
            fname = f'{calc_type}_{name}.xlsx'
            wb = openpyxl.Workbook()
            ws = wb.active
            ws.title = name

            if calc_type == 'RDF' and arr_y is not None:
                ws.append(['r (Å)', 'g(r)'])
                for r_val, g_val in zip(arr_x, arr_y):
                    ws.append([float(r_val), float(g_val)])
            elif calc_type == 'ADF':
                ws.append(['Angle (°)'])
                for a_val in arr_x:
                    ws.append([float(a_val)])
            else:
                ws.append(['Value'])
                for v in arr_x:
                    ws.append([float(v)])

            save_path = filedialog.asksaveasfilename(
                initialfile=fname,
                defaultextension='.xlsx',
                filetypes=[('Excel', '*.xlsx')])
            if save_path:
                wb.save(save_path)
                exported.append(save_path)

        if exported:
            messagebox.showinfo(self.t('conv_ok_title'),
                                self.t('export_ok').format('\n'.join(exported)))


# ═══════════════════════════════════════════════════════════════════════════════
# Dialogs
# ═══════════════════════════════════════════════════════════════════════════════
class _BaseDialog(tk.Toplevel):

    def __init__(self, parent, t_dict):
        super().__init__(parent)
        self._t = t_dict
        self.configure(bg=SURFACE)
        self.resizable(False, False)
        self.result = None
        self.grab_set()
        self.transient(parent)
        if IS_WIN:
            self.attributes('-toolwindow', True)

    def t(self, key):
        return self._t.get(key, key)

    def _row(self, parent, label, default='', w=11):
        row = tk.Frame(parent, bg=SURFACE)
        row.pack(fill='x', pady=5)
        tk.Label(row, text=label, bg=SURFACE, fg=FG2,
                 font=F_S, width=16, anchor='w').pack(side='left')
        var = tk.StringVar(value=default)
        ef = tk.Frame(row, bg=BORDER2, padx=1, pady=1)
        ef.pack(side='right')
        ttk.Entry(ef, textvariable=var, width=w, font=F_M).pack()
        return var


class _RdfParamsDialog(_BaseDialog):

    def __init__(self, parent, t_dict, elements):
        super().__init__(parent, t_dict)
        self.title(self.t('dlg_rdf_title'))

        outer = tk.Frame(self, bg=SURFACE)
        outer.pack(padx=28, pady=24)

        tk.Label(outer, text=self.t('dlg_rdf_title'),
                 bg=SURFACE, fg=FG, font=F_LB).pack(anchor='w')
        tk.Frame(outer, bg=BORDER, height=1).pack(fill='x', pady=(8, 12))

        elem_str = self.t('detected_elements').format('  ·  '.join(elements))
        tk.Label(outer, text=elem_str,
                 bg=SURFACE, fg=ACCENTL, font=F_S,
                 wraplength=280, justify='left').pack(anchor='w', pady=(0, 10))

        pairs = list(combinations_with_replacement(elements, 2))
        pairs_str = '  ·  '.join(f'{a}–{b}' for a, b in pairs)
        tk.Label(outer,
                 text=pairs_str,
                 bg=SURFACE, fg=FG2, font=F_XS,
                 wraplength=280, justify='left').pack(anchor='w', pady=(0, 10))

        tk.Frame(outer, bg=BORDER, height=1).pack(fill='x', pady=(0, 8))

        self.var_rmax = self._row(outer, self.t('rdf_rmax'), '10.0')
        self.var_dr   = self._row(outer, self.t('rdf_dr'),   '0.1')

        tk.Frame(outer, bg=BORDER, height=1).pack(fill='x', pady=(12, 10))
        btns = tk.Frame(outer, bg=SURFACE)
        btns.pack(fill='x')
        ttk.Button(btns, text=self.t('cancel'), style='Secondary.TButton',
                   command=self.destroy).pack(side='right', padx=(6, 0))
        ttk.Button(btns, text=self.t('ok'), style='Primary.TButton',
                   command=self._ok).pack(side='right')

    def _ok(self):
        try:
            rmax = float(self.var_rmax.get())
            dr   = float(self.var_dr.get())
            self.result = (dr, rmax)
            self.destroy()
        except ValueError as e:
            messagebox.showerror(self.t('error'), str(e), parent=self)


# ═══════════════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    app = BelloApp()
    app.mainloop()
