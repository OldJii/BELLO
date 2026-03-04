"""
BELLO GUI — Cross-platform tkinter interface (macOS / Windows / Linux)
with embedded matplotlib plots and Chinese / English language switching.
"""

import sys
import platform
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import threading
import os

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import BELLO_main
import Angle_Distribution_Function
import Radial_Pair_Distribution_Function
import Coordination_Heatmap
import vasp_converter

# ═══════════════════════════════════════════════════════════════════════
# Platform helpers
# ═══════════════════════════════════════════════════════════════════════
IS_MAC = platform.system() == 'Darwin'
IS_WIN = platform.system() == 'Windows'


def _pick_family():
    """Return the best available font family for the current platform."""
    if IS_MAC:
        return 'Helvetica Neue'
    if IS_WIN:
        return 'Segoe UI'
    return 'sans-serif'


_FF = _pick_family()
FONT_S     = (_FF, 9)
FONT_M     = (_FF, 10)
FONT_M_B   = (_FF, 10, 'bold')
FONT_L     = (_FF, 11, 'bold')
FONT_XL    = (_FF, 17, 'bold')
FONT_TITLE = (_FF, 20, 'bold')

# ═══════════════════════════════════════════════════════════════════════
# i18n — all UI strings
# ═══════════════════════════════════════════════════════════════════════
_STRINGS = {
    'en': {
        'app_title':        'BELLO — Bond Element Lattice Locality Order',
        'header':           'B.E.L.L.O',
        'subtitle':         'Local-order analysis for disordered systems',
        'lang_btn':         '中文',
        'sec_file':         'Input File',
        'browse':           'Browse…',
        'sec_vasp':         'VASP → XYZ Converter',
        'convert':          'Convert',
        'sec_params':       'Parameters',
        'auto_threshold':   'Automatic threshold',
        'threshold':        'Threshold',
        'tolerance':        'Tolerance',
        'divide_angle':     'Divide angle distributions (memory efficient)',
        'sec_cell':         'Unit Cell  (Å)',
        'sec_calc':         'Calculation',
        'mode_bello':       'BELLO Analysis',
        'mode_rdf':         'RDF  (Radial Distribution)',
        'mode_angle':       'Angle Distribution',
        'mode_coord':       'Coordination Heatmap',
        'btn_calc':         'Calculate',
        'ready':            'Ready',
        'running':          'Running…',
        'done':             'Done!',
        'error':            'Error',
        'plot_placeholder':  'Results will appear here after calculation.',
        'plot_tab':         'Plot',
        'dlg_rdf_title':    'RDF Parameters',
        'rdf_rmax':         'Maximum r :',
        'rdf_dr':           'Delta r :',
        'rdf_elem1':        'First element :',
        'rdf_elem2':        'Second element :',
        'ok':               'OK',
        'cancel':           'Cancel',
        'dlg_elem_title':   'Select Elements',
        'num_elements':     'Number of elements :',
        'element_n':        'Element {}:',
        'err_vasp':         'Please select a valid VASP file.',
        'err_file':         'Please select a valid {} file.',
        'err_float':        '{} must be a number, got: "{}"',
        'err_elem_empty':   'All element names are required.',
        'err_elem_name':    'Element names cannot be empty.',
        'conv_ok_title':    'Conversion Successful',
        'conv_ok_msg':      'Output: {}\nCell: {:.3f} × {:.3f} × {:.3f} Å\nFrames: {}',
    },
    'zh': {
        'app_title':        'BELLO — 键元素晶格局域有序度分析',
        'header':           'B.E.L.L.O',
        'subtitle':         '无序/非晶体系局域有序度分析工具',
        'lang_btn':         'English',
        'sec_file':         '输入文件',
        'browse':           '浏览…',
        'sec_vasp':         'VASP → XYZ 转换',
        'convert':          '转换',
        'sec_params':       '参数设置',
        'auto_threshold':   '自动确定阈值',
        'threshold':        '距离阈值',
        'tolerance':        '距离容差',
        'divide_angle':     '分段保存角度分布（节省内存）',
        'sec_cell':         '晶胞尺寸  (Å)',
        'sec_calc':         '计算',
        'mode_bello':       'BELLO 分析',
        'mode_rdf':         'RDF 径向分布函数',
        'mode_angle':       '角度分布函数',
        'mode_coord':       '配位热力图',
        'btn_calc':         '开始计算',
        'ready':            '就绪',
        'running':          '计算中…',
        'done':             '完成！',
        'error':            '错误',
        'plot_placeholder':  '计算完成后，结果图表将显示在此处。',
        'plot_tab':         '图表',
        'dlg_rdf_title':    'RDF 参数',
        'rdf_rmax':         '最大半径 r :',
        'rdf_dr':           '步长 Δr :',
        'rdf_elem1':        '第一种元素 :',
        'rdf_elem2':        '第二种元素 :',
        'ok':               '确定',
        'cancel':           '取消',
        'dlg_elem_title':   '选择元素',
        'num_elements':     '元素数量 :',
        'element_n':        '元素 {} :',
        'err_vasp':         '请选择有效的 VASP 文件。',
        'err_file':         '请选择有效的 {} 文件。',
        'err_float':        '{} 必须为数字，当前值: "{}"',
        'err_elem_empty':   '所有元素名称均为必填项。',
        'err_elem_name':    '元素名称不能为空。',
        'conv_ok_title':    '转换成功',
        'conv_ok_msg':      '输出: {}\n晶胞: {:.3f} × {:.3f} × {:.3f} Å\n帧数: {}',
    },
}

# ═══════════════════════════════════════════════════════════════════════
# Color palette  (Catppuccin Mocha inspired, tuned for readability)
# ═══════════════════════════════════════════════════════════════════════
C_BG        = '#1e1e2e'
C_SURFACE   = '#24243a'
C_CARD      = '#2a2a40'
C_CARD_HL   = '#313150'
C_FG        = '#cdd6f4'
C_FG_DIM    = '#7f849c'
C_ACCENT    = '#89b4fa'
C_ACCENT_HV = '#74c7ec'
C_BTN       = '#363654'
C_BTN_HV    = '#45456a'
C_ENTRY_BG  = '#313150'
C_BORDER    = '#45456a'
C_GREEN     = '#a6e3a1'
C_RED       = '#f38ba8'
C_YELLOW    = '#f9e2af'


# ═══════════════════════════════════════════════════════════════════════
# Main application
# ═══════════════════════════════════════════════════════════════════════
class BelloApp(tk.Tk):

    def __init__(self):
        super().__init__()

        self._lang = 'en'
        self._t = _STRINGS['en']

        self.title(self._t['app_title'])
        self.configure(bg=C_BG)
        self.minsize(1000, 680)
        self.geometry('1120x760')

        if IS_MAC:
            self.createcommand('tk::mac::Quit', self.destroy)

        self._apply_theme()
        self._build_ui()

    # ── helpers ──────────────────────────────────────────────────────
    def t(self, key):
        return self._t.get(key, key)

    def _switch_lang(self):
        self._lang = 'zh' if self._lang == 'en' else 'en'
        self._t = _STRINGS[self._lang]
        self.title(self.t('app_title'))
        self._rebuild_all()

    # ── theme ────────────────────────────────────────────────────────
    def _apply_theme(self):
        s = ttk.Style(self)
        s.theme_use('clam')

        s.configure('.', background=C_BG, foreground=C_FG, borderwidth=0,
                    font=FONT_M)
        s.configure('TFrame', background=C_BG)
        s.configure('Surface.TFrame', background=C_SURFACE)
        s.configure('Card.TFrame', background=C_CARD)

        s.configure('TLabel', background=C_BG, foreground=C_FG, font=FONT_M)
        s.configure('Title.TLabel', background=C_SURFACE, foreground=C_ACCENT,
                    font=FONT_TITLE)
        s.configure('Sub.TLabel', background=C_SURFACE, foreground=C_FG_DIM,
                    font=FONT_S)
        s.configure('Sec.TLabel', background=C_CARD, foreground=C_ACCENT,
                    font=FONT_L)
        s.configure('Card.TLabel', background=C_CARD, foreground=C_FG,
                    font=FONT_M)
        s.configure('Dim.TLabel', background=C_BG, foreground=C_FG_DIM,
                    font=FONT_M)
        s.configure('Status.TLabel', background=C_CARD, foreground=C_GREEN,
                    font=FONT_M_B)

        s.configure('TEntry', fieldbackground=C_ENTRY_BG, foreground=C_FG,
                    insertcolor=C_FG, borderwidth=1, relief='solid', font=FONT_M)

        for w in ('TCheckbutton', 'TRadiobutton'):
            s.configure(w, background=C_CARD, foreground=C_FG, font=FONT_M)
            s.map(w, background=[('active', C_CARD)])

        s.configure('Accent.TButton', background=C_ACCENT, foreground=C_BG,
                    font=FONT_L, padding=(18, 9))
        s.map('Accent.TButton',
              background=[('active', C_ACCENT_HV), ('pressed', C_ACCENT_HV)])

        s.configure('Tool.TButton', background=C_BTN, foreground=C_FG,
                    font=FONT_S, padding=(10, 5))
        s.map('Tool.TButton', background=[('active', C_BTN_HV)])

        s.configure('Lang.TButton', background=C_SURFACE, foreground=C_ACCENT,
                    font=FONT_S, padding=(8, 4))
        s.map('Lang.TButton', background=[('active', C_CARD)])

        s.configure('TNotebook', background=C_BG, borderwidth=0, tabmargins=(4, 4, 4, 0))
        s.configure('TNotebook.Tab', background=C_BTN, foreground=C_FG,
                    font=FONT_M, padding=(14, 6))
        s.map('TNotebook.Tab',
              background=[('selected', C_ACCENT)],
              foreground=[('selected', C_BG)])

        s.configure('Horizontal.TProgressbar', troughcolor=C_BTN,
                    background=C_ACCENT, borderwidth=0, thickness=6)

        s.configure('TSeparator', background=C_BORDER)

        s.configure('TCombobox', fieldbackground=C_ENTRY_BG, foreground=C_FG,
                    background=C_BTN, font=FONT_M)
        s.map('TCombobox', fieldbackground=[('readonly', C_ENTRY_BG)])

    # ── build / rebuild ──────────────────────────────────────────────
    def _rebuild_all(self):
        for w in self.winfo_children():
            w.destroy()
        self._build_ui()

    def _build_ui(self):
        self._build_header()

        body = ttk.Frame(self)
        body.pack(fill='both', expand=True, padx=14, pady=(0, 14))
        body.columnconfigure(1, weight=1)
        body.rowconfigure(0, weight=1)

        self._build_sidebar(body)
        self._build_plot_area(body)

    # ── header ───────────────────────────────────────────────────────
    def _build_header(self):
        hdr = ttk.Frame(self, style='Surface.TFrame')
        hdr.pack(fill='x', ipady=8)

        inner = ttk.Frame(hdr, style='Surface.TFrame')
        inner.pack(fill='x', padx=20)

        ttk.Label(inner, text=self.t('header'), style='Title.TLabel').pack(
            side='left')
        ttk.Label(inner, text=self.t('subtitle'), style='Sub.TLabel').pack(
            side='left', padx=(14, 0), pady=(8, 0))

        ttk.Button(inner, text=self.t('lang_btn'), style='Lang.TButton',
                   command=self._switch_lang).pack(side='right')

    # ── sidebar (left) ───────────────────────────────────────────────
    def _build_sidebar(self, parent):
        outer = ttk.Frame(parent)
        outer.grid(row=0, column=0, sticky='ns', padx=(0, 8))

        canvas = tk.Canvas(outer, bg=C_BG, highlightthickness=0, width=370)
        sb = ttk.Scrollbar(outer, orient='vertical', command=canvas.yview)
        self._scroll_frame = ttk.Frame(canvas)

        self._scroll_frame.bind(
            '<Configure>',
            lambda e: canvas.configure(scrollregion=canvas.bbox('all')))
        canvas.create_window((0, 0), window=self._scroll_frame, anchor='nw')
        canvas.configure(yscrollcommand=sb.set)

        canvas.pack(side='left', fill='both', expand=True)
        sb.pack(side='right', fill='y')

        self._sidebar_canvas = canvas
        self._bind_scroll(canvas)

        self._build_file_card(self._scroll_frame)
        self._build_vasp_card(self._scroll_frame)
        self._build_params_card(self._scroll_frame)
        self._build_cell_card(self._scroll_frame)
        self._build_calc_card(self._scroll_frame)

    def _bind_scroll(self, canvas):
        if IS_MAC:
            canvas.bind_all('<MouseWheel>',
                            lambda e: canvas.yview_scroll(-e.delta, 'units'))
        elif IS_WIN:
            canvas.bind_all('<MouseWheel>',
                            lambda e: canvas.yview_scroll(int(-e.delta / 120), 'units'))
        else:
            canvas.bind_all('<Button-4>',
                            lambda e: canvas.yview_scroll(-3, 'units'))
            canvas.bind_all('<Button-5>',
                            lambda e: canvas.yview_scroll(3, 'units'))

    # ── card factory ─────────────────────────────────────────────────
    def _card(self, parent, title_key):
        wrapper = ttk.Frame(parent, style='Card.TFrame')
        wrapper.pack(fill='x', padx=2, pady=5)
        body = ttk.Frame(wrapper, style='Card.TFrame')
        body.pack(fill='x', padx=16, pady=(10, 12))
        ttk.Label(body, text=self.t(title_key), style='Sec.TLabel').pack(
            anchor='w')
        ttk.Separator(body).pack(fill='x', pady=(4, 8))
        return body

    def _row(self, parent, label_text, default='', entry_width=18):
        r = ttk.Frame(parent, style='Card.TFrame')
        r.pack(fill='x', pady=3)
        ttk.Label(r, text=label_text, style='Card.TLabel').pack(side='left')
        var = tk.StringVar(value=default)
        ttk.Entry(r, textvariable=var, width=entry_width).pack(side='right')
        return var

    # ── file card ────────────────────────────────────────────────────
    def _build_file_card(self, parent):
        card = self._card(parent, 'sec_file')
        r = ttk.Frame(card, style='Card.TFrame')
        r.pack(fill='x', pady=2)
        self.url_var = tk.StringVar()
        ttk.Entry(r, textvariable=self.url_var, width=28).pack(
            side='left', fill='x', expand=True)
        ttk.Button(r, text=self.t('browse'), style='Tool.TButton',
                   command=self._browse_xyz).pack(side='right', padx=(6, 0))

    def _browse_xyz(self):
        p = filedialog.askopenfilename(
            title='XYZ',
            filetypes=[('XYZ', '*.xyz'), ('All files', '*')])
        if p:
            self.url_var.set(p)

    # ── vasp card ────────────────────────────────────────────────────
    def _build_vasp_card(self, parent):
        card = self._card(parent, 'sec_vasp')
        r = ttk.Frame(card, style='Card.TFrame')
        r.pack(fill='x', pady=2)
        self.vasp_var = tk.StringVar()
        ttk.Entry(r, textvariable=self.vasp_var, width=28).pack(
            side='left', fill='x', expand=True)
        ttk.Button(r, text=self.t('browse'), style='Tool.TButton',
                   command=self._browse_vasp).pack(side='right', padx=(6, 0))
        ttk.Button(card, text=self.t('convert'), style='Tool.TButton',
                   command=self._convert_vasp).pack(fill='x', pady=(6, 0))

    def _browse_vasp(self):
        p = filedialog.askopenfilename(
            title='VASP',
            filetypes=[('All files', '*'), ('VASP', '*.vasp')])
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

    # ── params card ──────────────────────────────────────────────────
    def _build_params_card(self, parent):
        card = self._card(parent, 'sec_params')
        self.auto_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(card, text=self.t('auto_threshold'),
                        variable=self.auto_var,
                        command=self._toggle_auto).pack(anchor='w', pady=(0, 4))
        self.trh_var = self._row(card, self.t('threshold'), '3.0')
        self.tlr_var = self._row(card, self.t('tolerance'), '0.5')
        self.sep_ang_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(card, text=self.t('divide_angle'),
                        variable=self.sep_ang_var).pack(anchor='w', pady=(6, 0))

    def _toggle_auto(self):
        st = 'disabled' if self.auto_var.get() else 'normal'
        for w in self._all_widgets():
            if isinstance(w, ttk.Entry):
                try:
                    vname = str(w.cget('textvariable'))
                    if vname in (str(self.trh_var), str(self.tlr_var)):
                        w.configure(state=st)
                except Exception:
                    pass

    def _all_widgets(self):
        stack = list(self.winfo_children())
        while stack:
            w = stack.pop()
            yield w
            stack.extend(w.winfo_children())

    # ── cell card ────────────────────────────────────────────────────
    def _build_cell_card(self, parent):
        card = self._card(parent, 'sec_cell')
        self.cell_x = self._row(card, 'X', '20.00000')
        self.cell_y = self._row(card, 'Y', '20.00000')
        self.cell_z = self._row(card, 'Z', '20.00000')

    # ── calc card ────────────────────────────────────────────────────
    def _build_calc_card(self, parent):
        card = self._card(parent, 'sec_calc')
        self.calc_mode = tk.StringVar(value='BELLO')
        modes = [
            ('mode_bello', 'BELLO'),
            ('mode_rdf',   'RDF'),
            ('mode_angle', 'ANGLE'),
            ('mode_coord', 'COORD'),
        ]
        for key, val in modes:
            ttk.Radiobutton(card, text=self.t(key), variable=self.calc_mode,
                            value=val).pack(anchor='w', pady=2)

        ttk.Button(card, text=self.t('btn_calc'), style='Accent.TButton',
                   command=self._on_calculate).pack(fill='x', pady=(14, 6))

        self.progress = ttk.Progressbar(card, mode='indeterminate',
                                        style='Horizontal.TProgressbar')
        self.progress.pack(fill='x', pady=(2, 4))

        self.status_var = tk.StringVar(value=self.t('ready'))
        ttk.Label(card, textvariable=self.status_var,
                  style='Status.TLabel').pack(anchor='w')

    # ── plot area (right) ────────────────────────────────────────────
    def _build_plot_area(self, parent):
        self.plot_notebook = ttk.Notebook(parent)
        self.plot_notebook.grid(row=0, column=1, sticky='nsew')

        ph = ttk.Frame(self.plot_notebook)
        self.plot_notebook.add(ph, text=f'  {self.t("plot_tab")}  ')
        ttk.Label(ph, text=self.t('plot_placeholder'),
                  style='Dim.TLabel',
                  font=FONT_L).place(relx=0.5, rely=0.5, anchor='center')

    def _show_figures(self, figs):
        for tab_id in self.plot_notebook.tabs():
            self.plot_notebook.forget(tab_id)
        for idx, fig in enumerate(figs):
            frame = ttk.Frame(self.plot_notebook)
            self.plot_notebook.add(
                frame, text=f'  {self.t("plot_tab")} {idx + 1}  ')
            canvas = FigureCanvasTkAgg(fig, master=frame)
            canvas.draw()
            tb_frame = ttk.Frame(frame)
            tb_frame.pack(side='bottom', fill='x')
            NavigationToolbar2Tk(canvas, tb_frame)
            canvas.get_tk_widget().pack(side='top', fill='both', expand=True)
        if figs:
            self.plot_notebook.select(0)

    # ── validation ───────────────────────────────────────────────────
    def _vfloat(self, value, name='value'):
        try:
            return float(value)
        except (ValueError, TypeError):
            raise ValueError(self.t('err_float').format(name, value))

    def _vfile(self, path, label='file'):
        if not path or not os.path.isfile(path):
            raise ValueError(self.t('err_file').format(label))

    # ── calculation dispatch ─────────────────────────────────────────
    def _on_calculate(self):
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
        self.progress.start(10)
        self.status_var.set(self.t('running'))

        def _worker():
            try:
                result = fn(*args)
                self.after(0, self._done, result)
            except Exception as exc:
                self.after(0, self._fail, exc)

        threading.Thread(target=_worker, daemon=True).start()

    def _done(self, figures):
        self.progress.stop()
        self.status_var.set(self.t('done'))
        if figures:
            self._show_figures(figures)

    def _fail(self, exc):
        self.progress.stop()
        self.status_var.set(self.t('error'))
        messagebox.showerror(self.t('error'), str(exc))

    def _pcb(self, stage, cur, total):
        pct = int(cur / max(total, 1) * 100)
        self.after(0, lambda: self.status_var.set(f'{stage}: {pct}%'))

    # ── BELLO ────────────────────────────────────────────────────────
    def _run_bello(self):
        url = self.url_var.get().strip()
        self._vfile(url, 'XYZ')
        cx = self._vfloat(self.cell_x.get(), 'X')
        cy = self._vfloat(self.cell_y.get(), 'Y')
        cz = self._vfloat(self.cell_z.get(), 'Z')
        auto = self.auto_var.get()
        trh = self.trh_var.get() if not auto else '2'
        tlr = self.tlr_var.get() if not auto else '1'
        self._threaded(BELLO_main.BELLO,
                       (url, cx, cy, cz, auto, trh, tlr,
                        self.sep_ang_var.get(), self._pcb))

    # ── RDF ──────────────────────────────────────────────────────────
    def _run_rdf(self):
        url = self.url_var.get().strip()
        self._vfile(url, 'XYZ')
        cx = self._vfloat(self.cell_x.get(), 'X')
        cy = self._vfloat(self.cell_y.get(), 'Y')
        cz = self._vfloat(self.cell_z.get(), 'Z')
        dlg = _RdfDialog(self, self._t)
        self.wait_window(dlg)
        if dlg.result is None:
            return
        dr, rmax, e1, e2 = dlg.result
        self._threaded(Radial_Pair_Distribution_Function.RDF,
                       (url, cx, cy, cz, dr, rmax, e1, e2, self._pcb))

    # ── Angle ────────────────────────────────────────────────────────
    def _run_angle(self):
        dlg = _ElemDialog(self, self._t)
        self.wait_window(dlg)
        if dlg.result is None:
            return
        self._threaded(Angle_Distribution_Function.sorter, (dlg.result,))

    # ── Coord ────────────────────────────────────────────────────────
    def _run_coord(self):
        dlg = _ElemDialog(self, self._t)
        self.wait_window(dlg)
        if dlg.result is None:
            return
        self._threaded(Coordination_Heatmap.coordination_heatmap, (dlg.result,))


# ═══════════════════════════════════════════════════════════════════════
# Dialogs
# ═══════════════════════════════════════════════════════════════════════
class _BaseDialog(tk.Toplevel):
    """Themed dialog base."""
    def __init__(self, parent, t_dict):
        super().__init__(parent)
        self._t = t_dict
        self.configure(bg=C_BG)
        self.resizable(False, False)
        self.result = None
        self.grab_set()
        self.transient(parent)
        if IS_WIN:
            self.attributes('-toolwindow', True)

    def t(self, key):
        return self._t.get(key, key)


class _RdfDialog(_BaseDialog):
    def __init__(self, parent, t_dict):
        super().__init__(parent, t_dict)
        self.title(self.t('dlg_rdf_title'))

        frm = ttk.Frame(self)
        frm.pack(padx=24, pady=20)

        fields = [
            ('rdf_rmax',  '10.0'),
            ('rdf_dr',    '0.1'),
            ('rdf_elem1', 'Ge'),
            ('rdf_elem2', 'Se'),
        ]
        self.vars = []
        for key, default in fields:
            row = ttk.Frame(frm)
            row.pack(fill='x', pady=4)
            ttk.Label(row, text=self.t(key), width=18).pack(side='left')
            v = tk.StringVar(value=default)
            ttk.Entry(row, textvariable=v, width=14).pack(side='right')
            self.vars.append(v)

        btns = ttk.Frame(frm)
        btns.pack(fill='x', pady=(14, 0))
        ttk.Button(btns, text=self.t('cancel'), style='Tool.TButton',
                   command=self.destroy).pack(side='right', padx=(6, 0))
        ttk.Button(btns, text=self.t('ok'), style='Accent.TButton',
                   command=self._ok).pack(side='right')

    def _ok(self):
        try:
            rmax = float(self.vars[0].get())
            dr = float(self.vars[1].get())
            e1 = self.vars[2].get().strip()
            e2 = self.vars[3].get().strip()
            if not e1 or not e2:
                raise ValueError(self.t('err_elem_name'))
            self.result = (dr, rmax, e1, e2)
            self.destroy()
        except ValueError as e:
            messagebox.showerror(self.t('error'), str(e), parent=self)


class _ElemDialog(_BaseDialog):
    def __init__(self, parent, t_dict):
        super().__init__(parent, t_dict)
        self.title(self.t('dlg_elem_title'))

        frm = ttk.Frame(self)
        frm.pack(padx=24, pady=20)

        row0 = ttk.Frame(frm)
        row0.pack(fill='x', pady=(0, 8))
        ttk.Label(row0, text=self.t('num_elements')).pack(side='left')
        self.num_var = tk.IntVar(value=2)
        cb = ttk.Combobox(row0, textvariable=self.num_var,
                          values=[2, 3, 4, 5, 6], width=5, state='readonly')
        cb.pack(side='right')
        cb.bind('<<ComboboxSelected>>', self._rebuild)

        self._fields = ttk.Frame(frm)
        self._fields.pack(fill='x')
        self.elem_vars = []
        self._rebuild()

        btns = ttk.Frame(frm)
        btns.pack(fill='x', pady=(14, 0))
        ttk.Button(btns, text=self.t('cancel'), style='Tool.TButton',
                   command=self.destroy).pack(side='right', padx=(6, 0))
        ttk.Button(btns, text=self.t('ok'), style='Accent.TButton',
                   command=self._ok).pack(side='right')

    def _rebuild(self, _event=None):
        for w in self._fields.winfo_children():
            w.destroy()
        self.elem_vars = []
        for i in range(self.num_var.get()):
            row = ttk.Frame(self._fields)
            row.pack(fill='x', pady=2)
            ttk.Label(row, text=self.t('element_n').format(i + 1),
                      width=12).pack(side='left')
            v = tk.StringVar()
            ttk.Entry(row, textvariable=v, width=10).pack(side='right')
            self.elem_vars.append(v)

    def _ok(self):
        elems = [v.get().strip() for v in self.elem_vars]
        if any(not e for e in elems):
            messagebox.showerror(self.t('error'), self.t('err_elem_empty'),
                                 parent=self)
            return
        self.result = elems
        self.destroy()


# ═══════════════════════════════════════════════════════════════════════
# Entry point
# ═══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    app = BelloApp()
    app.mainloop()
