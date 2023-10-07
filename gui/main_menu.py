import os

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import QUrl, Qt
from PyQt6.QtGui import QPixmap, QPalette
from PyQt6.QtWebEngineWidgets import QWebEngineView
from PyQt6.QtWidgets import QFileDialog, QTableWidgetItem, QAbstractItemView, QLineEdit
import sys


class Ui_Group(object):
    def getOutputPath(self):
        file = open('quantifySpliceVariants.nf', 'r')
        final = file.readlines()
        final = final[-3].split("out = ")[1]
        final = final.split("\n")[0]
        self.loadCsvOnOpen(final + "/hallo.tsv")
        self.pixmap = QPixmap(final + "/Bananas.svg")
        self.label.setPixmap(self.pixmap)
        self.webView.setUrl(QUrl("file:///" + final + "/result.pdf"))
        file.close()

    def open_file_dialog(self):
        file_dialog = QFileDialog()
        filename, _ = file_dialog.getOpenFileName()
        if filename:
            self.webView.setUrl(QUrl("file:///" + filename.replace('\\', '/')))

    def openImage(self):
        fname = QFileDialog.getOpenFileName(Group, "Open File", os.getenv('HOME'),
                                            "SVG Files (*.svg);;PNG Files (*.png);;Jpg Files (*.jpg)")
        # Open The Image
        if fname:
            self.pixmap = QPixmap(fname[0])
            # Add Pic to label
            self.label.setPixmap(self.pixmap)

    def loadCsvOnOpen(self, fileName):
        if fileName:
            if ".csv" in fileName:
                df = pd.read_csv(fileName)
                header = df.iloc[0]
                df = df[1:]
                self.tableView.setColumnCount(len(df.columns))
                self.tableView.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                    for j in range(len(df.columns)):
                        m = QTableWidgetItem(header[j])
                        self.tableView.setHorizontalHeaderItem(j, m)
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
                header = df.iloc[0]
                self.tableView.setColumnCount(len(df.columns))
                self.tableView.setRowCount(len(df.index))
                for i in range(len(df.index)):
                    for j in range(len(df.columns)):
                        self.tableView.setItem(i, j, QTableWidgetItem(str(df.iat[i, j])))
                        for j in range(len(df.columns)):
                            m = QTableWidgetItem(header[j])
                            self.tableView.setHorizontalHeaderItem(j, m)
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
        fileName, _ = QFileDialog.getOpenFileName(Group, "Open CSV",
                                                  os.getenv('Home'), "CSV Files (*.csv);; TSV Files (*.tsv);; "
                                                                     "TXT Files (*.txt) ;; TAB Files (*.tab)")
        if fileName:
            self.loadCsvOnOpen(fileName)

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
    def setupUi(self, Group):
        self.window = QtWidgets.QMainWindow()
        Group.setObjectName("Group")
        Group.showFullScreen()
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

        self.deTab = QtWidgets.QTabWidget()
        self.deTab.setGeometry(QtCore.QRect(0, 0, 791, 631))
        self.deTab.setObjectName("tabWidget")
        self.de = QtWidgets.QWidget()
        self.de.setObjectName("de")
        self.horizontalLayoutWidget3 = QtWidgets.QWidget(parent=self.de)
        self.horizontalLayoutWidget3.resize(1500, 900)
        self.horizontalLayoutWidget3.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout3 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget3)
        self.horizontalLayout3.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout3.setObjectName("horizontalLayout")
        self.verticalLayout3 = QtWidgets.QVBoxLayout()
        self.query2 = QLineEdit()
        self.query2.setPlaceholderText("Search...")
        self.query2.textChanged.connect(self.search)
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
        self.label3 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget3)
        self.label3.setObjectName("label")
        self.pixmap3 = QPixmap("Logo.png")
        self.label3.setPixmap(self.pixmap)
        self.horizontalLayout3.addLayout(self.verticalLayout3)
        self.horizontalLayout3.addWidget(self.label3)
        self.tabWidget.addTab(self.de, "")

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
        self.actionProject = QtGui.QAction(parent=Group)
        self.actionProject.setObjectName("actionProject")
        self.actionClose = QtGui.QAction(parent=Group)
        self.actionClose.setObjectName("actionClose")
        self.actionClose.triggered.connect(quit)
        self.isoformAllAction = QtGui.QAction(parent=Group)
        self.isoformAllAction.setObjectName("isoformAllAction")
        self.menuOpen.addAction(self.actionProject)
        self.menuOpen.addSeparator()
        self.menuOpen.addAction(self.actionCSV)
        self.menuOpen.addAction(self.actionStatistics)
        self.menuOpen.addAction(self.actionReport)
        self.menuFile.addAction(self.menuOpen.menuAction())
        self.menuFile.addAction(self.actionClose)
        self.menuIsoform.addAction(self.isoformAllAction)
        self.menuIsoform.addSeparator()
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuIsoform.menuAction())

        self.getOutputPath()

        font.setPointSize(10)

        self.retranslateUi(Group)
        QtCore.QMetaObject.connectSlotsByName(Group)

    def retranslateUi(self, Group):
        _translate = QtCore.QCoreApplication.translate
        Group.setWindowTitle(_translate("Group", "Main Menu"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuOpen.setTitle(_translate("MainWindow", "Open"))
        self.actionCSV.setText(_translate("MainWindow", "CSV"))
        self.actionReport.setText(_translate("MainWindow", "Report"))
        self.actionStatistics.setText(_translate("MainWindow", "Statistics"))
        self.actionProject.setText(_translate("MainWindow", "Project"))
        self.actionClose.setText(_translate("MainWindow", "Close"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.statistics), _translate("MainWindow", "Statistics"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.report), _translate("MainWindow", "Report"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.de), _translate("MainWindow", "Differential Gene Expression"))


if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    Group = QtWidgets.QMainWindow()
    ui = Ui_Group()
    ui.setupUi(Group)
    Group.show()
    sys.exit(app.exec())