%\LoadClass{article}
\NeedsTeXFormat{LaTeX2e}
\LoadClass[UTF8,adobefonts,a4paper,fntef]{ctexart}
\RequirePackage{amsmath,amsfonts,amssymb,bm}%数学相关宏包支持
\RequirePackage{graphicx}%提供图形支持
\RequirePackage[left=2cm,right=2cm,top=2cm,bottom=2.5cm]{geometry}%页面设置宏包
\RequirePackage{titlesec}%标题格式支持
\RequirePackage{fancyhdr}
\RequirePackage{enumitem}
\RequirePackage{enumerate}%列表环境支持
\RequirePackage{listings}%插入代码
\RequirePackage{url}%插入超链接
\RequirePackage{xcolor}%颜色
\RequirePackage{fontspec}%中文字体
\RequirePackage{float}%图形位置
\RequirePackage{subfigure}%提供图形支持
\RequirePackage{siunitx,array,booktabs}
\RequirePackage{csvsimple}%可使用csv文件数据


\linespread{1.2}%设置行间距
\setcounter{secnumdepth}{2}%
%设置各级标题格式
\CTEXsetup[name={实验,}]{part}
\CTEXsetup[format={\zihao{3}\bf\kaishu\raggedright}]{part}
\CTEXsetup[name={(,)}]{section}
\CTEXsetup[number={\chinese{section}}]{section}
\CTEXsetup[format={\zihao{4}\bf\kaishu\raggedright}]{section}
\CTEXsetup[number={},titleformat={\zihao{-4}\fangsong}]{subsection}
\newlist{clause}{enumerate}{1}
\setlist[clause]{fullwidth,itemindent=2em,label=(\arabic*)}

\setmainfont{Times New Roman}%设置主英文字体为 Times New Roman
\setCJKmainfont[ItalicFont={Adobe Fangsong Std R}]{Adobe Song Std} %中文主字体对应宋体
\setmonofont{Consolas}%等宽字体（代码）

%代码段样式设置
\lstset{
  language=[90]Fortran,
  extendedchars=false,
  frame=tb,
  aboveskip=3mm,
  belowskip=3mm,
  xleftmargin=2em,
  xrightmargin=2em,
  showstringspaces=false,
  columns=flexible,
  framerule=1pt,
  rulecolor=\color{gray!35},
  backgroundcolor=\color{gray!5},
  basicstyle={\footnotesize\ttfamily},
  numbers=left,
  numberstyle=\color{magenta},
  keywordstyle=\bfseries,%代码关键字的颜色为粗体
  commentstyle=\color[HTML]{898B8E},%浅灰色的注释
  stringstyle=\color{red},% 代码字符串为红色
  identifierstyle=\color[HTML]{10069F},
  breaklines=true,
  breakatwhitespace=true,
  tabsize=4,
}

\newcommand{\sid}{161xx1yyy}% 学号
\newcommand{\name}{Yvonne Wong}% 姓名
\newcommand{\dat}{\today}
\newcommand{\numm}{1}

%封面生成
\newcommand{\makecover}
  {\thispagestyle{empty}
  \begin{titlepage}\vspace{3em}
    \includegraphics[width=0.14\textwidth]{njulogo.pdf}\vspace{2.6em}
    \begin{center}
      \includegraphics[width=0.5\textwidth]{njuname.pdf}\\[2.8em]
      \zihao{1}\textbf{第{\numm}次实习报告}
      \vfill
      \begin{minipage}{0.8\textwidth}
        \begin{description}
        \zihao{3}\fangsong
        \item[\textbf{学 ~~~号}]      \uline{\hfill\textit{\sid}\hfill~}
        \item[\textbf{学生姓名}]      \uline{\hfill\textit{\name}\hfill~}
        \item[\textbf{提交日期}]      \uline{\hfill\textit{\dat}\hfill~}
        \end{description}
      \end{minipage}
      \vfill
    \end{center}
  \end{titlepage}\cleardoublepage}
