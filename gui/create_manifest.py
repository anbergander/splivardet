import csv
import os

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from parameter import Ui_Dialog


class Ui_createManifest(object):
    def openWindow(self, path):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Dialog()
        self.ui.setupUi(self.window)
        self.ui.plainTextEdit_6.setPlainText(path + ".tsv")
        self.window.show()
    def updateTable(self):
        self.spinBox.setDisabled(True)
        self.pushButton_3.setDisabled(True)
        self.pushButton.setDisabled(False)
        self.tableView.setDisabled(False)
        self.tableView.setColumnCount(4)
        col_headers = ['Sample', 'Group', 'Batch', 'Path']
        self.tableView.setHorizontalHeaderLabels(col_headers)
        self.tableView.setRowCount(self.spinBox.value())

    def save_sheet(self):
        for row in range(self.tableView.rowCount()):
            for column in range(self.tableView.columnCount()):
                item = self.tableView.item(row, column)
                if item is not None:
                    if column == 0:
                        if "barcode" in item.text():
                            continue
                        else:
                            QMessageBox.information(self.centralwidget, "Error", "Sample name must contain barcode + a number")
                            return
                    if column == 3:
                        isdir = os.path.isfile(item.text())
                        if isdir:
                            continue
                        else:
                            QMessageBox.information(self.centralwidget, "Error", "Must be an existing path")
                            return
                else:
                    QMessageBox.information(self.centralwidget, "Error", "Table must be filled completely")
                    return
        path = QFileDialog.getSaveFileName(self.centralwidget, 'Save TSV', os.getenv('HOME'), 'TSV (*.tsv)')
        if path[0] != '':
            with open(path[0], 'w') as csv_file:
                writer = csv.writer(csv_file, delimiter="\t", lineterminator='\n')
                for row in range(self.tableView.rowCount()):
                    row_data = []
                    for column in range(self.tableView.columnCount()):
                        item = self.tableView.item(row, column)
                        if item is not None:
                            row_data.append(item.text())
                        else:
                            row_data.append('')
                    writer.writerow(row_data)
        self.openWindow(path[0])

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(687, 594)
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pushButton = QtWidgets.QPushButton(parent=self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(540, 530, 131, 41))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setDisabled(True)
        self.pushButton.clicked.connect(self.save_sheet)
        self.pushButton.clicked.connect(MainWindow.close)
        self.pushButton.setObjectName("pushButton")
        font = QtGui.QFont()
        font.setPointSize(10)
        self.pushButton.setFont(font)
        self.pushButton.setObjectName("pushButton")
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("next-step-28.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.pushButton.setIcon(icon)

        self.verticalLayoutWidget = QtWidgets.QWidget(parent=self.centralwidget)
        self.verticalLayoutWidget.setGeometry(QtCore.QRect(39, 10, 611, 511))
        self.verticalLayoutWidget.setObjectName("verticalLayoutWidget")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.verticalLayoutWidget)
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_3.addWidget(self.label_2)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_3)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(parent=self.verticalLayoutWidget)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.spinBox = QtWidgets.QSpinBox(parent=self.verticalLayoutWidget)
        self.spinBox.setMinimum(2)
        self.spinBox.setObjectName("spinBox")
        self.horizontalLayout.addWidget(self.spinBox)
        self.pushButton_3 = QtWidgets.QPushButton(parent=self.verticalLayoutWidget)
        self.pushButton_3.setObjectName("pushButton_3")
        icon2 = QtGui.QIcon()
        icon2.addPixmap(QtGui.QPixmap("../../ressources/add.png"), QtGui.QIcon.Mode.Normal,
                       QtGui.QIcon.State.Off)
        self.pushButton_3.setIcon(icon2)
        self.pushButton_3.clicked.connect(self.updateTable)

        self.horizontalLayout.addWidget(self.pushButton_3)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.tableView = QtWidgets.QTableWidget(parent=self.verticalLayoutWidget)
        self.tableView.setObjectName("tableView")
        self.tableView.setDisabled(True)
        self.verticalLayout.addWidget(self.tableView)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 687, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Create Manifest"))
        self.pushButton.setText(_translate("MainWindow", "Continue"))
        self.label_2.setText(_translate("MainWindow", "Group"))
        self.label.setText(_translate("MainWindow", "Number of Samples"))
        self.pushButton_3.setText(_translate("MainWindow", "Add"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_createManifest()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec())