from PyQt6 import QtCore, QtGui, QtWidgets
from create_manifest import Ui_createManifest



class Ui_LoadManifest(object):
    def openWindowLoad(self):
        from load_mainfest import Ui_MainWindow
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.window)
        self.window.show()
    def openWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_createManifest()
        self.ui.setupUi(self.window)
        self.window.show()

    def setupUi(self, LoadDialog):
        self.window = self.window = QtWidgets.QMainWindow()
        LoadDialog.setObjectName("loadReport")
        LoadDialog.setFixedSize(400, 215)
        font = QtGui.QFont()
        font.setFamily("Calibri")
        LoadDialog.setFont(font)

        self.label = QtWidgets.QLabel(parent=LoadDialog)
        self.label.setGeometry(QtCore.QRect(20, 10, 361, 91))
        font.setPointSize(11)
        self.label.setFont(font)
        self.label.setObjectName("label")


        self.centralwidget = QtWidgets.QWidget(LoadDialog)
        self.centralwidget.setObjectName("centralwidget")
        LoadDialog.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(LoadDialog)


        self.pushButton = QtWidgets.QPushButton(self.centralwidget, clicked = lambda: self.openWindow())
        self.pushButton.clicked.connect(LoadDialog.close)
        self.pushButton.setGeometry(QtCore.QRect(250, 150, 121, 41))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap("next-step-28.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.pushButton.setIcon(icon)
        font.setPointSize(10)
        self.pushButton.setFont(font)

        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget, clicked=lambda: self.openWindowLoad())
        self.pushButton_2.clicked.connect(LoadDialog.close)
        self.pushButton_2.setGeometry(QtCore.QRect(20, 150, 121, 41))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("load_39552.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.pushButton_2.setIcon(icon1)
        self.pushButton_2.setFont(font)

        self.menubar.setGeometry(QtCore.QRect(0, 0, 515, 21))
        self.menubar.setObjectName("menubar")
        LoadDialog.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(LoadDialog)
        self.statusbar.setObjectName("statusbar")
        LoadDialog.setStatusBar(self.statusbar)

        self.retranslateUi(LoadDialog)
        QtCore.QMetaObject.connectSlotsByName(LoadDialog)

    def retranslateUi(self, LoadDialog):
        _translate = QtCore.QCoreApplication.translate
        self.pushButton.setText(_translate("LoadDialog", "Continue"))
        self.pushButton_2.setText(_translate("LoadDialog", "Load"))
        LoadDialog.setWindowTitle(_translate("LoadDialog", "Manifest"))
        self.label.setText(_translate("LoadDialog",
                                      "<html><head/><body><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">Load Manifest or</span></p><p align=\"center\"><span style=\" font-size:14pt; font-weight:600;\">create a new Manifest?</span></p></body></html>"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    window = QtWidgets.QMainWindow()
    ui = Ui_LoadManifest()
    ui.setupUi(window)
    window.show()
    sys.exit(app.exec())