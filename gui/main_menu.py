import os

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QUrl, Qt
from PyQt6.QtGui import QPixmap, QPalette
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QFileDialog, QTableWidgetItem, QAbstractItemView, QLineEdit
import sys

from gui.de_dialog import Ui_DE
from gui.go_dialog import Ui_GO
from gui.report_de import Ui_ReportDe


class Ui_Group(object):
    def openCreateReport(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_ReportDe()
        self.ui.setupUi(self.window)
        self.window.show()
    def openWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_GO()
        self.ui.setupUi(self.window)
        self.window.show()
    def openWindowDe(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_DE()
        self.ui.setupUi(self.window)
        self.window.show()
    def getOutputPath(self):
        file = open('quantifySpliceVariants.nf', 'r')
        final = file.readlines()
        final = final[-3].split("out = ")[1]
        final = final.split("\n")[0]
        self.pixmap = QPixmap(final + "/figures/voilin.svg")
        self.label.setPixmap(self.pixmap)
        self.webView.setUrl(QUrl("file:///" + final + "/result.pdf"))
        file.close()

    def open_file_dialog(self):
        file_dialog = QFileDialog()
        filename, _ = file_dialog.getOpenFileName()
        if filename:
            self.webView.setUrl(QUrl("file:///" + filename.replace('\\', '/')))

    def openImage(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        if fname:
            self.pixmap = QPixmap(fname[0])
            self.label.setPixmap(self.pixmap)

    def openHeatmap(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap3 = QPixmap(fname[0])
            # Add Pic to label
            self.label_4.setPixmap(self.pixmap3)
    def openPCA(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap2 = QPixmap(fname[0])
            # Add Pic to label
            self.label_2.setPixmap(self.pixmap2)
    def openVolcano(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap4 = QPixmap(fname[0])
            # Add Pic to label
            self.label_3.setPixmap(self.pixmap4)
    def openBarplot(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap5 = QPixmap(fname[0])
            # Add Pic to label
            self.label_5.setPixmap(self.pixmap5)

    def openLinePlot(self):
        fname = QFileDialog.getOpenFileName(self.centralwidget, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap6 = QPixmap(fname[0])
            # Add Pic to label
            self.label_6.setPixmap(self.pixmap6)

    def loadCsvOnOpen(self, fileName):
        if fileName:
            if ".csv" in fileName:
                df = pd.read_csv(fileName)
                df = df[1:]
                self.tableView.setColumnCount(len(df.columns))
                self.tableView.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView.selectRow(0)
                self.tableView.resizeColumnsToContents()
                self.tableView.resizeRowsToContents()
                col_headers = ['Gene', 'Isoform', 'Chromosome', 'p-Value', 'GO terms']
                self.tableView.setHorizontalHeaderLabels(col_headers)
                self.tableView.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView.setColumnWidth(0, 100)
                self.tableView.setColumnWidth(1, 100)
                self.tableView.setColumnWidth(2, 100)
                self.tableView.setColumnWidth(3, 100)
                self.tableView.setColumnWidth(4, 100)
            else:
                df = pd.read_table(fileName, sep='\t', lineterminator="\n")
                self.tableView.setColumnCount(len(df.columns))
                self.tableView.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView.selectRow(0)
                self.tableView.resizeColumnsToContents()
                self.tableView.resizeRowsToContents()
                col_headers = ['Gene', 'Isoform', 'Chromosome', 'p-Value', 'GO terms']
                self.tableView.setHorizontalHeaderLabels(col_headers)
                self.tableView.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView.setColumnWidth(0, 100)
                self.tableView.setColumnWidth(1, 100)
                self.tableView.setColumnWidth(2, 100)
                self.tableView.setColumnWidth(3, 100)
                self.tableView.setColumnWidth(4, 100)

    def loadCsv(self):
        fileName, _ = QFileDialog.getOpenFileName(self.centralwidget, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        if fileName:
            self.loadCsvOnOpen(fileName)
    def loadCsvDe(self):
        fileName, _ = QFileDialog.getOpenFileName(self.centralwidget, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        if fileName:
            self.loadDeTable(fileName)

    def loadCsvGO(self):
        fileName, _ = QFileDialog.getOpenFileName(self.centralwidget, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        if fileName:
            self.loadGoTable(fileName)


    def loadGoTable(self, fileName):
        if fileName:
            if ".csv" in fileName:
                df = pd.read_csv(fileName)
                df = df[1:]
                self.tableView3.setColumnCount(len(df.columns))
                self.tableView3.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView3.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView3.selectRow(0)
                self.tableView3.resizeColumnsToContents()
                self.tableView3.resizeRowsToContents()
                col_headers = ['Term', 'FDR', 'ES', 'NES']
                self.tableView3.setHorizontalHeaderLabels(col_headers)
                self.tableView3.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView3.setColumnWidth(0, 100)
                self.tableView3.setColumnWidth(1, 100)
                self.tableView3.setColumnWidth(2, 100)
                self.tableView3.setColumnWidth(3, 100)
                self.tableView3.setColumnWidth(4, 100)
            else:
                df = pd.read_table(fileName, sep='\t', lineterminator="\n")
                self.tableView3.setColumnCount(len(df.columns))
                self.tableView3.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView3.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView3.selectRow(0)
                self.tableView3.resizeColumnsToContents()
                self.tableView3.resizeRowsToContents()
                col_headers = ['Term', 'FDR', 'ES', 'NES']
                self.tableView3.setHorizontalHeaderLabels(col_headers)
                self.tableView3.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView3.setColumnWidth(0, 100)
                self.tableView3.setColumnWidth(1, 100)
                self.tableView3.setColumnWidth(2, 100)
                self.tableView3.setColumnWidth(3, 100)
                self.tableView3.setColumnWidth(4, 100)

    def loadDeTable(self, fileName):
        if fileName:
            if ".csv" in fileName:
                df = pd.read_csv(fileName)
                df = df[1:]
                self.tableView2.setColumnCount(len(df.columns))
                self.tableView2.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView2.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView2.selectRow(0)
                self.tableView2.resizeColumnsToContents()
                self.tableView2.resizeRowsToContents()
                col_headers = ['Gene', 'Base Mean', 'log2-Fold Change', 'lfcSE', 'Stat', 'P-value', 'Padj']
                self.tableView2.setHorizontalHeaderLabels(col_headers)
                self.tableView2.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView2.setColumnWidth(0, 100)
                self.tableView2.setColumnWidth(1, 100)
                self.tableView2.setColumnWidth(2, 100)
                self.tableView2.setColumnWidth(3, 100)
                self.tableView2.setColumnWidth(4, 100)
            else:
                df = pd.read_table(fileName, sep='\t', lineterminator="\n")
                self.tableView2.setColumnCount(len(df.columns))
                self.tableView2.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView2.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                self.tableView2.selectRow(0)
                self.tableView2.resizeColumnsToContents()
                self.tableView2.resizeRowsToContents()
                col_headers = ['Gene', 'Base Mean', 'log2-Fold Change', 'lfcSE', 'Stat', 'P-value', 'Padj']
                self.tableView2.setHorizontalHeaderLabels(col_headers)
                self.tableView2.sortItems(3, order=Qt.SortOrder.AscendingOrder)
                self.tableView2.setColumnWidth(0, 100)
                self.tableView2.setColumnWidth(1, 100)
                self.tableView2.setColumnWidth(2, 100)
                self.tableView2.setColumnWidth(3, 100)
                self.tableView2.setColumnWidth(4, 100)

    def search(self, s):
        # Clear current selection.
        self.tableView.setCurrentItem(None)

        if not s:
            # Empty string, don't search.
            return

        matching_items = self.tableView.findItems(s, Qt.MatchFlag.MatchContains)
        if matching_items:
            # We have found something.
            item = matching_items[0]  # Take the first.
            self.tableView.setCurrentItem(item)

    def search2(self, s):
        # Clear current selection.
        self.tableView2.setCurrentItem(None)

        if not s:
            # Empty string, don't search.
            return

        matching_items = self.tableView2.findItems(s, Qt.MatchFlag.MatchContains)
        if matching_items:
            # We have found something.
            item = matching_items[0]  # Take the first.
            self.tableView2.setCurrentItem(item)

    def search3(self, s):
        # Clear current selection.
        self.tableView3.setCurrentItem(None)

        if not s:
            # Empty string, don't search.
            return

        matching_items = self.tableView3.findItems(s, Qt.MatchFlag.MatchContains)
        if matching_items:
            # We have found something.
            item = matching_items[0]  # Take the first.
            self.tableView3.setCurrentItem(item)

    def setupUi(self, Group):
        self.window = QtWidgets.QMainWindow()
        Group.setObjectName("Group")
        Group.resize(900, 900)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        font.setBold(True)
        font.setWeight(75)
        Group.setFont(font)

        self.centralwidget = QtWidgets.QWidget(Group)
        self.centralwidget.setObjectName("centralwidget")
        Group.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(Group)

        self.layout = QtWidgets.QVBoxLayout()

        self.tabWidget = QtWidgets.QTabWidget()
        self.tabWidget.setGeometry(QtCore.QRect(0, 0, 791, 631))
        self.tabWidget.setObjectName("tabWidget")
        self.statistics = QtWidgets.QWidget()
        self.statistics.setObjectName("statistics")
        self.horizontalLayoutWidget = QtWidgets.QWidget(parent=self.statistics)
        self.horizontalLayoutWidget.resize(1500, 900)
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.query = QLineEdit()
        self.query.setPlaceholderText("Search...")
        self.query.textChanged.connect(self.search)
        self.root = os.path.dirname(sys.argv[0])
        self.delimit = '\t'
        Group.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.fileName = ""
        self.tableView = QtWidgets.QTableWidget(parent=self.horizontalLayoutWidget)
        self.tableView.setObjectName("tableView")
        self.tableView.setCornerButtonEnabled(False)
        self.tableView.setShowGrid(True)
        self.tableView.horizontalHeader().setBackgroundRole(QPalette.ColorRole.Window)
        self.tableView.setDropIndicatorShown(True)
        self.tableView.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.loadCsvOnOpen("examplemanifest.tsv")
        self.verticalLayout.addWidget(self.query)
        self.verticalLayout.addWidget(self.tableView)
        self.label = QtWidgets.QLabel(parent=self.horizontalLayoutWidget)
        self.label.setObjectName("label")
        self.pixmap = QPixmap("Logo.png")
        self.label.setPixmap(self.pixmap)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.horizontalLayout.addWidget(self.label)
        self.tabWidget.addTab(self.statistics, "")

        self.de = QtWidgets.QWidget()
        self.de.setObjectName("de")
        self.horizontalLayoutWidget3 = QtWidgets.QWidget(parent=self.de)
        self.horizontalLayoutWidget3.resize(2000, 900)
        self.horizontalLayoutWidget3.setObjectName("horizontalLayoutWidget")
        self.verticalLayout3 = QtWidgets.QVBoxLayout(self.horizontalLayoutWidget3)
        self.query2 = QLineEdit()
        self.query2.setPlaceholderText("Search...")
        self.query2.textChanged.connect(self.search2)
        self.root2 = os.path.dirname(sys.argv[0])
        self.delimit2 = '\t'
        Group.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.fileName2 = ""
        self.tableView2 = QtWidgets.QTableWidget(parent=self.horizontalLayoutWidget3)
        self.tableView2.setObjectName("tableView")
        self.tableView2.setCornerButtonEnabled(False)
        self.tableView2.setShowGrid(True)
        self.tableView2.horizontalHeader().setBackgroundRole(QPalette.ColorRole.Window)
        self.tableView2.setDropIndicatorShown(True)
        self.tableView2.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.verticalLayout3.addWidget(self.query2)
        self.verticalLayout3.addWidget(self.tableView2)
        self.tabWidget.addTab(self.de, "")

        self.GeneExpression = QtWidgets.QWidget()
        self.GeneExpression.setObjectName("GeneExpression")
        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(parent=self.GeneExpression)
        self.horizontalLayoutWidget_2.resize(2000, 850)
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_4 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_2)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_2.addWidget(self.label_4)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_2)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.label_3 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_2)
        self.label_3.setObjectName("label_3")
        self.verticalLayout_2.addWidget(self.label_3)
        self.horizontalLayout_2.addLayout(self.verticalLayout_2)
        self.tabWidget.addTab(self.GeneExpression, "")

        self.goTable = QtWidgets.QWidget()
        self.goTable.setObjectName("de")
        self.horizontalLayoutWidget_4 = QtWidgets.QWidget(parent=self.goTable)
        self.horizontalLayoutWidget_4.resize(2000, 900)
        self.horizontalLayoutWidget_4.setObjectName("horizontalLayoutWidget")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.horizontalLayoutWidget_4)
        self.query3 = QLineEdit()
        self.query3.setPlaceholderText("Search...")
        self.query3.textChanged.connect(self.search3)
        self.root3 = os.path.dirname(sys.argv[0])
        self.delimit3 = '\t'
        Group.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.fileName3 = ""
        self.tableView3 = QtWidgets.QTableWidget(parent=self.horizontalLayoutWidget_4)
        self.tableView3.setObjectName("tableView")
        self.tableView3.setCornerButtonEnabled(False)
        self.tableView3.setShowGrid(True)
        self.tableView3.horizontalHeader().setBackgroundRole(QPalette.ColorRole.Window)
        self.tableView3.setDropIndicatorShown(True)
        self.tableView3.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.verticalLayout_3.addWidget(self.query3)
        self.verticalLayout_3.addWidget(self.tableView3)
        self.tabWidget.addTab(self.goTable, "")

        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(parent=self.tab)
        self.horizontalLayoutWidget_3.resize(1900, 800)
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_5 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_3)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_3.addWidget(self.label_5)
        self.label_6 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_3)
        self.label_6.setObjectName("label_6")
        self.horizontalLayout_3.addWidget(self.label_6)
        self.tabWidget.addTab(self.tab, "")

        self.report = QtWidgets.QWidget()
        self.report.setObjectName("report")
        self.verticalLayoutWidget = QtWidgets.QWidget(parent=self.report)
        self.verticalLayoutWidget.resize(2000, 900)
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.layoutReport = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.layoutReport.setObjectName("layoutReport")
        self.webView = QWebEngineView(parent=self.verticalLayoutWidget)
        self.webView.settings().setAttribute(self.webView.settings().WebAttribute.PluginsEnabled, True)
        self.webView.settings().setAttribute(self.webView.settings().WebAttribute.PdfViewerEnabled, True)
        self.layoutReport.addWidget(self.webView)
        self.tabWidget.addTab(self.report, "")

        self.layout.addWidget(self.tabWidget)

        self.centralwidget.setLayout(self.layout)

        Group.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=Group)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 848, 21))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(parent=self.menubar)
        self.menuFile.setObjectName("menuResult")
        self.menuOpen = QtWidgets.QMenu(parent=self.menuFile)
        self.menuIsoform = QtWidgets.QMenu(parent=self.menubar)
        self.menuIsoform.setObjectName("menuIsoform")
        self.menuDE = QtWidgets.QMenu(parent=self.menubar)
        self.menuDE.setObjectName("menuDE")
        self.menuPerform = QtWidgets.QMenu(parent=self.menuDE)
        Group.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=Group)
        self.statusbar.setObjectName("statusbar")
        Group.setStatusBar(self.statusbar)
        self.actionCSV = QtGui.QAction(parent=Group)
        self.actionCSV.setObjectName("actionCSV")
        self.actionCSV.triggered.connect(self.loadCsv)
        self.actionReport = QtGui.QAction(parent=Group)
        self.actionReport.setObjectName("actionReport")
        self.actionReport.triggered.connect(self.open_file_dialog)
        self.actionStatistics = QtGui.QAction(parent=Group)
        self.actionStatistics.setObjectName("actionStatistics")
        self.actionStatistics.triggered.connect(self.openImage)
        self.actionDeTable = QtGui.QAction(parent=Group)
        self.actionDeTable.setObjectName("actionDeTable")
        self.actionDeTable.triggered.connect(self.loadCsvDe)
        self.actionDeHeatmap = QtGui.QAction(parent=Group)
        self.actionDeHeatmap.setObjectName("actionDeHeatmap")
        self.actionDeHeatmap.triggered.connect(self.openHeatmap)
        self.actionDePCA = QtGui.QAction(parent=Group)
        self.actionDePCA.setObjectName("actionDePCA")
        self.actionDePCA.triggered.connect(self.openPCA)
        self.actionDeVolcano = QtGui.QAction(parent=Group)
        self.actionDeVolcano.setObjectName("actionDePCA")
        self.actionDeVolcano.triggered.connect(self.openVolcano)
        self.actionGOTable = QtGui.QAction(parent=Group)
        self.actionGOTable.setObjectName("actionDePCA")
        self.actionGOTable.triggered.connect(self.loadCsvGO)
        self.actionGOgraph = QtGui.QAction(parent=Group)
        self.actionGOgraph.setObjectName("actionDePCA")
        self.actionGOgraph.triggered.connect(self.openLinePlot)
        self.actionGOgraphbar = QtGui.QAction(parent=Group)
        self.actionGOgraphbar.setObjectName("actionDePCA")
        self.actionGOgraphbar.triggered.connect(self.openBarplot)
        self.actionDeReport = QtGui.QAction(parent=Group)
        self.actionDeReport.setObjectName("actionDeReport")
        self.actionDeReport.triggered.connect(self.openCreateReport)
        self.actionProject = QtGui.QAction(parent=Group)
        self.actionProject.setObjectName("actionProject")
        self.actionClose = QtGui.QAction(parent=Group)
        self.actionClose.setObjectName("actionClose")
        self.actionClose.triggered.connect(Group.close)
        self.actionDe = QtGui.QAction(parent=Group)
        self.actionDe.setObjectName("actionDe")
        self.actionDe.triggered.connect(self.openWindowDe)
        self.actionGo = QtGui.QAction(parent=Group)
        self.actionGo.setObjectName("actionGo")
        self.actionGo.triggered.connect(self.openWindow)
        self.isoformAllAction = QtGui.QAction(parent=Group)
        self.isoformAllAction.setObjectName("isoformAllAction")
        #self.menuOpen.addAction(self.actionProject)
        #self.menuOpen.addSeparator()
        self.menuOpen.addAction(self.actionCSV)
        self.menuOpen.addAction(self.actionStatistics)
        self.menuOpen.addSeparator()
        self.menuOpen.addAction(self.actionDeTable)
        self.menuOpen.addAction(self.actionDeHeatmap)
        self.menuOpen.addAction(self.actionDePCA)
        self.menuOpen.addAction(self.actionDeVolcano)
        self.menuOpen.addSeparator()
        self.menuOpen.addAction(self.actionGOTable)
        self.menuOpen.addAction(self.actionGOgraphbar)
        self.menuOpen.addAction(self.actionGOgraph)
        self.menuOpen.addSeparator()
        self.menuOpen.addAction(self.actionReport)
        self.menuFile.addAction(self.menuOpen.menuAction())
        self.menuFile.addAction(self.actionClose)
        self.menuIsoform.addAction(self.isoformAllAction)
        self.menuIsoform.addSeparator()
        self.menuPerform.addAction(self.actionDe)
        self.menuPerform.addAction(self.actionGo)
        self.menuDE.addAction(self.menuPerform.menuAction())
        self.menuDE.addSeparator()
        self.menuDE.addAction(self.actionDeReport)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuIsoform.menuAction())
        self.menubar.addAction(self.menuDE.menuAction())

        self.getOutputPath()

        font.setPointSize(10)

        self.retranslateUi(Group)
        QtCore.QMetaObject.connectSlotsByName(Group)

    def retranslateUi(self, Group):
        _translate = QtCore.QCoreApplication.translate
        Group.setWindowTitle(_translate("Group", "Main Menu"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuOpen.setTitle(_translate("MainWindow", "Open"))
        self.menuIsoform.setTitle(_translate("MainWindow", "Isoform"))
        self.menuDE.setTitle(_translate("MainWindow", "Differential Gene Expression"))
        self.menuPerform.setTitle(_translate("MainWindow", "Perform"))
        self.actionCSV.setText(_translate("MainWindow", "Isoform Table"))
        self.actionReport.setText(_translate("MainWindow", "Report"))
        self.actionStatistics.setText(_translate("MainWindow", "Isoform Plots"))
        self.actionDeTable.setText(_translate("MainWindow", "Differential Gene Expression Table"))
        self.actionDeHeatmap.setText(_translate("MainWindow", "Differential Gene Expression Heatmap"))
        self.actionProject.setText(_translate("MainWindow", "Project"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.actionDe.setText(_translate("MainWindow", "Differential Expression Analysis"))
        self.actionGo.setText(_translate("MainWindow", "GO Terms"))
        self.actionDePCA.setText(_translate("MainWindow", "Differential Gene Expression PCA"))
        self.actionDeVolcano.setText(_translate("MainWindow", "Differential Gene Expression Volcano"))
        self.actionGOgraphbar.setText(_translate("MainWindow", "GO Terms Barplot"))
        self.actionGOgraph.setText(_translate("MainWindow", "GO Terms Line Plot"))
        self.actionGOTable.setText(_translate("MainWindow", "GO Terms Table"))
        self.actionDeReport.setText(_translate("MainWindow", "Create Report for DE Analysis"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.statistics), _translate("MainWindow", "Isoform"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.report), _translate("MainWindow", "Report"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.de), _translate("MainWindow", "Differential Gene Expression Table"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.GeneExpression), _translate("MainWindow", "Differential Gene Expression"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab),_translate("MainWindow", "GO Terms"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.goTable), _translate("MainWindow", "GO Terms Table"))


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_Group()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())