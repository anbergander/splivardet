import os

from PyQt6 import QtCore, QtGui, QtWidgets
import multiprocessing

from PyQt6.QtWidgets import QFileDialog


class Ui_Dialog(object):
    def writeConfig(self):
        file = open('../quantifySpliceVariants.nf', 'w')
        s = ["params { \n", "\n", "threads = " + str(self.spinBox.value()) + " \n", "desaltIndex = " + self.plainTextEdit.toPlainText()
             + "\n", "genomeFasta = " + self.plainTextEdit_8.toPlainText() + " \n", "annotationsGtf = " + self.plainTextEdit_5.toPlainText() + " \n",
             "indexDir = " + self.plainTextEdit_11.toPlainText() + " \n", "dirToFastqs = " + self.plainTextEdit_4.toPlainText() + " \n",
             "flairManifest = " + self.plainTextEdit_6.toPlainText() + " \n", "out = " + self.plainTextEdit_7.toPlainText() + " \n", " \n", "}"]
        file.writelines(s)
        file.close()

    def browseGtf(self):
        fname = QFileDialog.getOpenFileName(self.horizontalLayoutWidget_4, "Open file", os.getenv('Home'), 'Gtf (*.gtf)')
        self.plainTextEdit_5.setPlainText(fname[0])

    def browseIndex(self):
        fname = QFileDialog.getExistingDirectory(self.horizontalLayoutWidget_3, "Select Directory", os.getenv('Home'))
        self.plainTextEdit.setPlainText(fname)
    def browseFasta(self):
        fname = QFileDialog.getOpenFileName(self.horizontalLayoutWidget, "Open file", os.getenv('Home'), 'Fasta (*.fasta) ;; Fasta (*.fa) ;; Genomefasta (*.genomefasta)')
        self.plainTextEdit_8.setPlainText(fname[0])

    def browseDirectory(self):
        fname = QFileDialog.getExistingDirectory(self.horizontalLayoutWidget_3, "Select Directory", os.getenv('Home'))
        self.plainTextEdit_4.setPlainText(fname)

    def browseDirectoryForOutput(self):
        fname = QFileDialog.getExistingDirectory(self.horizontalLayoutWidget_6, "Select Directory", os.getenv('Home'))
        self.plainTextEdit_7.setPlainText(fname)


    def browseDirectoryForIndex(self):
        fname = QFileDialog.getExistingDirectory(self.horizontalLayoutWidget_6, "Select Directory", os.getenv('Home'))
        self.plainTextEdit_11.setPlainText(fname)

    def checkIfFilled(self):
        if (self.plainTextEdit.toPlainText() == '' or self.plainTextEdit_4.toPlainText() == '' or self.plainTextEdit_5.toPlainText() == ''
                or self.plainTextEdit_6.toPlainText() == '' or self.plainTextEdit_8.toPlainText() == '' or self.plainTextEdit_11.toPlainText() == ''):
            self.pushButton_9.setDisabled(True)
        else:
            self.pushButton_9.setDisabled(False)

    def openWindow(self):
        self.writeConfig()
        from start_run import Ui_runDialog
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_runDialog()
        self.ui.setupUi(self.window)
        self.window.show()

    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.setFixedSize(545, 730)

        self.checkBox = QtWidgets.QCheckBox(parent=Dialog)
        self.checkBox.setGeometry(QtCore.QRect(35, 140, 150, 25))
        self.checkBox.setObjectName("checkBox")
        self.horizontalLayoutWidget = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget.setGeometry(QtCore.QRect(40, 90, 461, 41))
        self.horizontalLayoutWidget.setObjectName("horizontalLayoutWidget")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget)
        self.horizontalLayout.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(parent=self.horizontalLayoutWidget)
        self.label.setObjectName("label")
        self.label.setMinimumWidth(125)
        self.horizontalLayout.addWidget(self.label)
        self.plainTextEdit = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget)
        self.plainTextEdit.setObjectName("plainTextEdit")
        self.plainTextEdit.setMinimumWidth(225)
        self.plainTextEdit.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit.setDisabled(True)
        self.horizontalLayout.addWidget(self.plainTextEdit)
        self.pushButton = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.clicked.connect(self.browseIndex)
        self.pushButton.setMaximumWidth(100)
        self.horizontalLayout.addWidget(self.pushButton)

        self.horizontalLayoutWidget_2 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_2.setGeometry(QtCore.QRect(40, 180, 461, 41))
        self.horizontalLayoutWidget_2.setObjectName("horizontalLayoutWidget_2")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_2)
        self.horizontalLayout_2.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_2)
        self.label_2.setObjectName("label_2")
        self.label_2.setMaximumWidth(125)
        self.horizontalLayout_2.addWidget(self.label_2)
        self.spinBox = QtWidgets.QSpinBox(parent=self.horizontalLayoutWidget_2)
        self.spinBox.setObjectName("spinBox")
        self.spinBox.setMinimum(1)
        self.spinBox.setValue(10)
        self.spinBox.setMinimumWidth(250)
        self.spinBox.setMinimumHeight(40)
        self.spinBox.setMaximum(multiprocessing.cpu_count()-1)
        self.horizontalLayout_2.addWidget(self.spinBox)

        self.horizontalLayoutWidget_3 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_3.setGeometry(QtCore.QRect(40, 250, 461, 41))
        self.horizontalLayoutWidget_3.setObjectName("horizontalLayoutWidget_3")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_3)
        self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.label_4 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_3)
        self.label_4.setObjectName("label_4")
        self.label_4.setMinimumWidth(125)
        self.horizontalLayout_4.addWidget(self.label_4)
        self.plainTextEdit_4 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_3)
        self.plainTextEdit_4.setObjectName("plainTextEdit_4")
        self.plainTextEdit_4.setMinimumWidth(225)
        self.plainTextEdit_4.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_4.setDisabled(True)
        self.horizontalLayout_4.addWidget(self.plainTextEdit_4)
        self.pushButton_4 = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget_3)
        self.pushButton_4.setObjectName("pushButton_4")
        self.pushButton_4.setMaximumWidth(100)
        self.pushButton_4.clicked.connect(self.browseDirectory)
        self.horizontalLayout_4.addWidget(self.pushButton_4)

        self.horizontalLayoutWidget_4 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_4.setGeometry(QtCore.QRect(40, 320, 461, 41))
        self.horizontalLayoutWidget_4.setObjectName("horizontalLayoutWidget_4")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_4)
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.label_5 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_4)
        self.label_5.setObjectName("label_5")
        self.label_5.setMinimumWidth(125)
        self.horizontalLayout_5.addWidget(self.label_5)
        self.plainTextEdit_5 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_4)
        self.plainTextEdit_5.setObjectName("plainTextEdit_5")
        self.plainTextEdit_5.setMinimumWidth(225)
        self.plainTextEdit_5.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_5.setDisabled(True)
        self.horizontalLayout_5.addWidget(self.plainTextEdit_5)
        self.pushButton_5 = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget_4)
        self.pushButton_5.setObjectName("pushButton_5")
        self.pushButton_5.setMaximumWidth(100)
        self.pushButton_5.clicked.connect(self.browseGtf)
        self.horizontalLayout_5.addWidget(self.pushButton_5)

        self.horizontalLayoutWidget_5 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_5.setGeometry(QtCore.QRect(40, 390, 461, 41))
        self.horizontalLayoutWidget_5.setObjectName("horizontalLayoutWidget_5")
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_5)
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_6 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_5)
        self.label_6.setObjectName("label_6")
        self.label_6.setMinimumWidth(125)
        self.horizontalLayout_6.addWidget(self.label_6)
        self.plainTextEdit_6 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_5)
        self.plainTextEdit_6.setObjectName("plainTextEdit_6")
        self.plainTextEdit_6.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_6.setMinimumWidth(225)
        self.plainTextEdit_6.setDisabled(True)
        self.plainTextEdit_6.setDisabled(True)
        self.horizontalLayout_6.addWidget(self.plainTextEdit_6)

        self.horizontalLayoutWidget_6 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_6.setGeometry(QtCore.QRect(40, 460, 461, 41))
        self.horizontalLayoutWidget_6.setObjectName("horizontalLayoutWidget_6")
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_6)
        self.horizontalLayout_7.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_7 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_6)
        self.label_7.setObjectName("label_7")
        self.label_7.setMinimumWidth(125)
        self.horizontalLayout_7.addWidget(self.label_7)
        self.plainTextEdit_7 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_6)
        self.plainTextEdit_7.setObjectName("plainTextEdit_7")
        self.plainTextEdit_7.setMinimumWidth(225)
        self.plainTextEdit_7.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_7.setDisabled(True)
        self.horizontalLayout_7.addWidget(self.plainTextEdit_7)
        self.pushButton_7 = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget_6)
        self.pushButton_7.setObjectName("pushButton_7")
        self.pushButton_7.setMaximumWidth(100)
        self.pushButton_7.clicked.connect(self.browseDirectoryForOutput)
        self.horizontalLayout_7.addWidget(self.pushButton_7)

        self.horizontalLayoutWidget_7 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_7.setGeometry(QtCore.QRect(40, 530, 461, 41))
        self.horizontalLayoutWidget_7.setObjectName("horizontalLayoutWidget_6")
        self.horizontalLayout_8 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_7)
        self.horizontalLayout_8.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_8.setObjectName("horizontalLayout_7")
        self.label_9 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_7)
        self.label_9.setObjectName("label_7")
        self.label_9.setMinimumWidth(125)
        self.horizontalLayout_8.addWidget(self.label_9)
        self.plainTextEdit_8 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_7)
        self.plainTextEdit_8.setObjectName("plainTextEdit_7")
        self.plainTextEdit_8.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_8.setMinimumWidth(225)
        self.plainTextEdit_8.setDisabled(True)
        self.horizontalLayout_8.addWidget(self.plainTextEdit_8)
        self.pushButton_8 = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget_7)
        self.pushButton_8.setObjectName("pushButton_7")
        self.pushButton_8.setMaximumWidth(100)
        self.pushButton_8.clicked.connect(self.browseFasta)
        self.horizontalLayout_8.addWidget(self.pushButton_8)

        self.horizontalLayoutWidget_8 = QtWidgets.QWidget(parent=Dialog)
        self.horizontalLayoutWidget_8.setGeometry(QtCore.QRect(40, 600, 461, 41))
        self.horizontalLayoutWidget_8.setObjectName("horizontalLayoutWidget_6")
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout(self.horizontalLayoutWidget_8)
        self.horizontalLayout_9.setContentsMargins(0, 0, 0, 0)
        self.horizontalLayout_9.setObjectName("horizontalLayout_7")
        self.label_10 = QtWidgets.QLabel(parent=self.horizontalLayoutWidget_8)
        self.label_10.setObjectName("label_7")
        self.label_10.setMinimumWidth(125)
        self.horizontalLayout_9.addWidget(self.label_10)
        self.plainTextEdit_11 = QtWidgets.QPlainTextEdit(parent=self.horizontalLayoutWidget_8)
        self.plainTextEdit_11.setObjectName("plainTextEdit_11")
        self.plainTextEdit_11.textChanged.connect(self.checkIfFilled)
        self.plainTextEdit_11.setMinimumWidth(225)
        self.plainTextEdit_11.setDisabled(True)
        self.horizontalLayout_9.addWidget(self.plainTextEdit_11)
        self.pushButton_10 = QtWidgets.QPushButton(parent=self.horizontalLayoutWidget_8)
        self.pushButton_10.setObjectName("pushButton_7")
        self.pushButton_10.setMaximumWidth(100)
        self.pushButton_10.clicked.connect(self.browseDirectoryForIndex)
        self.horizontalLayout_9.addWidget(self.pushButton_10)


        self.label_8 = QtWidgets.QLabel(parent=Dialog)
        self.label_8.setGeometry(QtCore.QRect(210, 30, 231, 20))
        self.label_8.setObjectName("label_8")
        font = QtGui.QFont()
        font.setFamily("Calibri")
        font.setPointSize(16)
        self.label_8.setFont(font)
        self.pushButton_9 = QtWidgets.QPushButton(parent=Dialog)
        self.pushButton_9.setGeometry(QtCore.QRect(400, 670, 121, 41))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap("next-step-28.png"), QtGui.QIcon.Mode.Normal, QtGui.QIcon.State.Off)
        self.pushButton_9.setIcon(icon1)
        self.pushButton_9.setDisabled(True)
        self.pushButton_9.clicked.connect(self.openWindow)
        self.pushButton_9.clicked.connect(Dialog.close)
        self.pushButton_9.setObjectName("pushButton_9")

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Parameter for Config"))
        self.checkBox.setText(_translate("Dialog", "Build DeSALT index"))
        self.label.setText(_translate("Dialog", "Reference Genome"))
        self.pushButton.setText(_translate("Dialog", "Browse"))
        self.plainTextEdit.setPlaceholderText(_translate("Dialog", "Path to the reference genome/ index"))
        self.label_2.setText(_translate("Dialog", "Threads"))
        self.label_4.setText(_translate("Dialog", "FastQ Files"))
        self.plainTextEdit_4.setPlaceholderText(_translate("Dialog", "Path to FASTQ Files"))
        self.pushButton_4.setText(_translate("Dialog", "Browse"))
        self.label_5.setText(_translate("Dialog", "gtf File"))
        self.plainTextEdit_5.setPlaceholderText(_translate("Dialog", "Path to gtf File for annotation"))
        self.pushButton_5.setText(_translate("Dialog", "Browse"))
        self.label_6.setText(_translate("Dialog", "Flair Manifest"))
        self.plainTextEdit_6.setPlaceholderText(_translate("Dialog", "Path to Flair Manifest"))
        self.label_7.setText(_translate("Dialog", "Output"))
        self.plainTextEdit_7.setPlaceholderText(_translate("Dialog", "Optional"))
        self.pushButton_7.setText(_translate("Dialog", "Browse"))
        self.plainTextEdit_8.setPlaceholderText(_translate("Dialog", "Path to Genome Fasta"))
        self.label_8.setText(_translate("Dialog", "Select Parameters"))
        self.label_9.setText(_translate("Dialog", "Genome Fasta"))
        self.pushButton_8.setText(_translate("Dialog", "Browse"))
        self.pushButton_10.setText(_translate("Dialog", "Browse"))
        self.pushButton_9.setText(_translate("Dialog", "Continue"))
        self.plainTextEdit_11.setPlaceholderText(_translate("Dialog", "Path to Index"))
        self.label_10.setText(_translate("Dialog", "Index Directory"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QMainWindow()
    ui = Ui_Dialog()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec())