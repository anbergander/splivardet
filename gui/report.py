import os
from datetime import datetime

from PyQt6 import QtCore, QtGui, QtWidgets
from PyQt6.QtCore import pyqtSignal, QObject, QRunnable, pyqtSlot, QThreadPool, QUrl
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from pdfrw import PdfReader
from pdfrw.buildxobj import pagexobj
from pdfrw.toreportlab import makerl
from reportlab.pdfgen.canvas import Canvas

from main_menu import Ui_Group


class WorkerSignals(QObject):
    error = pyqtSignal(str)
    file_saved_as = pyqtSignal(str)

class PdfGenerator(QRunnable):
    def __init__(self, data, folder):
        super().__init__()
        self.folder = folder
        self.data = data
        self.signals = WorkerSignals()

    @pyqtSlot()
    def run(self):
        try:
            outfile = self.folder + "//" + "result.pdf"

            template = PdfReader("report.pdf", decompress=False).pages[0]
            template_obj = pagexobj(template)

            canvas = Canvas(outfile)

            xobj_name = makerl(canvas, template_obj)
            canvas.doForm(xobj_name)

            ystart = 720

            today = datetime.today()
            canvas.drawString(450, ystart, today.strftime('%F'))

            canvas.drawString(180, ystart - 45, self.data['experiment'])
            canvas.drawString(180, ystart - 67, self.data['reference'])
            canvas.drawString(180, ystart - 89, self.data['gtf'])


            canvas.drawString(180, ystart - 163, self.data['trimming'])
            canvas.drawString(180, ystart - 185, self.data['trimmingParam'])

            canvas.drawString(180, ystart - 247, self.data['dataprocessing'])
            canvas.drawString(180, ystart - 270, self.data['dataprocessingParam'])

            canvas.drawString(180, ystart - 336, self.data['alignment'])
            canvas.drawString(180, ystart - 357, self.data['alignmentParam'])

            canvas.drawString(180, ystart - 423, self.data['postprocessing'])
            canvas.drawString(180, ystart - 445, self.data['postprocessingParam'])

            canvas.drawString(180, ystart - 508, self.data['isoform'])
            canvas.drawString(180, ystart - 531, self.data['isoformParam'])

            canvas.drawString(180, ystart - 595, self.data['statistics'])
            canvas.drawString(180, ystart - 617, self.data['statisticsParam'])

            canvas.save()

        except Exception as e:
            self.signals.error.emit(str(e))
            return

        self.signals.file_saved_as.emit(outfile)

class Ui_Report(object):
    def openWindow(self):
        self.window = QtWidgets.QMainWindow()
        self.ui = Ui_Group()
        self.ui.setupUi(self.window)
        self.ui.webView.setUrl(QUrl("file:///" + self.path + "//" + "result.pdf"))
        self.window.show()
    def generate(self):
        fname = QFileDialog.getExistingDirectory(self.centralwidget, "Select Directory", os.getenv('Home'))
        fname = fname.replace('\\', '/')
        self.path = fname
        if fname:
            self.pushButton.setDisabled(True)
            data = {
                'experiment': self.textEdit_6.toPlainText(),
                'gtf': self.textEdit_2.toPlainText(),
                'reference': self.textEdit.toPlainText(),
                'trimming': self.textEdit_7.toPlainText(),
                'trimmingParam': self.textEdit_8.toPlainText(),
                'dataprocessing': self.textEdit_17.toPlainText(),
                'dataprocessingParam': self.textEdit_18.toPlainText(),
                'alignment': self.textEdit_19.toPlainText(),
                'alignmentParam': self.textEdit_20.toPlainText(),
                'postprocessing': self.textEdit_9.toPlainText(),
                'postprocessingParam': self.textEdit_10.toPlainText(),
                'isoform': self.textEdit_13.toPlainText(),
                'isoformParam': self.textEdit_14.toPlainText(),
                'statistics': self.textEdit_15.toPlainText(),
                'statisticsParam': self.textEdit_16.toPlainText(),
            }
            g = PdfGenerator(data, fname)
            g.signals.file_saved_as.connect(self.generated)
            g.signals.error.connect(print)
            self.threadpool.start(g)
            self.openWindow()

    def generated(self, outfile):
        Dialog.close()
        self.pushButton.setDisabled(False)
        try:
            os.startfile(outfile)
            Dialog.close()
        except Exception:
            QMessageBox.information(Dialog, "Finished", "PDF has been generated")

    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.setFixedSize(600, 780)
        self.path = ""
        self.threadpool = QThreadPool()
        self.centralwidget = QtWidgets.QWidget(parent=MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.verticalLayoutWidget_4 = QtWidgets.QWidget(parent=self.centralwidget)
        self.verticalLayoutWidget_4.setGeometry(QtCore.QRect(25, 15, 550, 730))
        self.verticalLayoutWidget_4.setObjectName("verticalLayoutWidget_4")
        self.verticalLayout_14 = QtWidgets.QVBoxLayout(self.verticalLayoutWidget_4)
        self.verticalLayout_14.setContentsMargins(0, 0, 0, 0)
        self.verticalLayout_14.setObjectName("verticalLayout_14")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem)
        self.label_23 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        font = QtGui.QFont()
        font.setPointSize(16)
        font.setFamily("Calibri")
        font.setBold(True)
        font.setWeight(75)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.horizontalLayout_3.addWidget(self.label_23)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Policy.Expanding, QtWidgets.QSizePolicy.Policy.Minimum)
        self.horizontalLayout_3.addItem(spacerItem1)
        self.verticalLayout_14.addLayout(self.horizontalLayout_3)
        self.horizontalLayout_11 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_11.setObjectName("horizontalLayout_11")
        self.label_6 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_6.setObjectName("label_6")
        self.label_6.setMinimumWidth(150)
        self.horizontalLayout_11.addWidget(self.label_6)
        self.textEdit_6 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_6.setObjectName("textEdit_6")
        self.textEdit_6.setMaximumHeight(30)
        self.horizontalLayout_11.addWidget(self.textEdit_6)
        self.verticalLayout_14.addLayout(self.horizontalLayout_11)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label.setObjectName("label")
        self.label.setMinimumWidth(150)
        self.horizontalLayout.addWidget(self.label)
        self.textEdit = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit.setObjectName("textEdit")
        self.textEdit.setMaximumHeight(30)
        self.horizontalLayout.addWidget(self.textEdit)
        self.verticalLayout_14.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_2.setObjectName("label_2")
        self.label_2.setMinimumWidth(150)
        self.horizontalLayout_2.addWidget(self.label_2)
        self.textEdit_2 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_2.setObjectName("textEdit_2")
        self.textEdit_2.setMaximumHeight(30)
        self.horizontalLayout_2.addWidget(self.textEdit_2)
        self.verticalLayout_14.addLayout(self.horizontalLayout_2)
        self.horizontalLayout_6 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_6.setObjectName("horizontalLayout_6")
        self.label_8 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_8.setObjectName("label_8")
        self.label_8.setMinimumWidth(150)
        self.horizontalLayout_6.addWidget(self.label_8)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_9 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_9.setObjectName("label_9")
        self.label_9.setMinimumWidth(150)
        self.verticalLayout_3.addWidget(self.label_9)
        self.label_10 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_10.setObjectName("label_10")
        self.label_10.setMinimumWidth(150)
        self.verticalLayout_3.addWidget(self.label_10)
        self.horizontalLayout_6.addLayout(self.verticalLayout_3)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout()
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.textEdit_7 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_7.setObjectName("textEdit_7")
        self.textEdit_7.setMinimumHeight(30)
        self.textEdit_7.setDisabled(True)
        self.verticalLayout_4.addWidget(self.textEdit_7)
        self.textEdit_8 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_8.setObjectName("textEdit_8")
        self.textEdit_8.setMinimumHeight(30)
        self.textEdit_8.setDisabled(True)
        self.verticalLayout_4.addWidget(self.textEdit_8)
        self.horizontalLayout_6.addLayout(self.verticalLayout_4)
        self.verticalLayout_14.addLayout(self.horizontalLayout_6)
        self.horizontalLayout_12 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_12.setObjectName("horizontalLayout_12")
        self.label_24 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_24.setObjectName("label_24")
        self.label_24.setMinimumWidth(150)
        self.horizontalLayout_12.addWidget(self.label_24)
        self.verticalLayout_13 = QtWidgets.QVBoxLayout()
        self.verticalLayout_13.setObjectName("verticalLayout_13")
        self.label_25 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_25.setObjectName("label_25")
        self.label_25.setMinimumWidth(150)
        self.verticalLayout_13.addWidget(self.label_25)
        self.label_26 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_26.setObjectName("label_26")
        self.label_26.setMinimumWidth(150)
        self.verticalLayout_13.addWidget(self.label_26)
        self.horizontalLayout_12.addLayout(self.verticalLayout_13)
        self.verticalLayout_15 = QtWidgets.QVBoxLayout()
        self.verticalLayout_15.setObjectName("verticalLayout_15")
        self.textEdit_17 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_17.setObjectName("textEdit_17")
        self.textEdit_17.setMinimumHeight(30)
        self.textEdit_17.setDisabled(True)
        self.verticalLayout_15.addWidget(self.textEdit_17)
        self.textEdit_18 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_18.setObjectName("textEdit_18")
        self.textEdit_18.setDisabled(True)
        self.textEdit_18.setMinimumHeight(30)
        self.verticalLayout_15.addWidget(self.textEdit_18)
        self.horizontalLayout_12.addLayout(self.verticalLayout_15)
        self.verticalLayout_14.addLayout(self.horizontalLayout_12)
        self.horizontalLayout_13 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_13.setObjectName("horizontalLayout_13")
        self.label_27 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_27.setObjectName("label_27")
        self.label_27.setMinimumWidth(150)
        self.horizontalLayout_13.addWidget(self.label_27)
        self.verticalLayout_16 = QtWidgets.QVBoxLayout()
        self.verticalLayout_16.setObjectName("verticalLayout_16")
        self.label_28 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_28.setObjectName("label_28")
        self.label_28.setMinimumWidth(150)
        self.verticalLayout_16.addWidget(self.label_28)
        self.label_29 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_29.setObjectName("label_29")
        self.label_29.setMinimumWidth(150)
        self.verticalLayout_16.addWidget(self.label_29)
        self.horizontalLayout_13.addLayout(self.verticalLayout_16)
        self.verticalLayout_17 = QtWidgets.QVBoxLayout()
        self.verticalLayout_17.setObjectName("verticalLayout_17")
        self.textEdit_19 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_19.setObjectName("textEdit_19")
        self.textEdit_19.setMinimumHeight(30)
        self.textEdit_19.setDisabled(True)
        self.verticalLayout_17.addWidget(self.textEdit_19)
        self.textEdit_20 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_20.setObjectName("textEdit_20")
        self.textEdit_20.setMinimumHeight(30)
        self.textEdit_20.setDisabled(True)
        self.verticalLayout_17.addWidget(self.textEdit_20)
        self.horizontalLayout_13.addLayout(self.verticalLayout_17)
        self.verticalLayout_14.addLayout(self.horizontalLayout_13)
        self.horizontalLayout_7 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_7.setObjectName("horizontalLayout_7")
        self.label_11 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_11.setObjectName("label_11")
        self.label_11.setMinimumWidth(150)
        self.horizontalLayout_7.addWidget(self.label_11)
        self.verticalLayout_5 = QtWidgets.QVBoxLayout()
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.label_12 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_12.setObjectName("label_12")
        self.label_12.setMinimumWidth(150)
        self.verticalLayout_5.addWidget(self.label_12)
        self.label_13 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_13.setObjectName("label_13")
        self.label_13.setMinimumWidth(150)
        self.verticalLayout_5.addWidget(self.label_13)
        self.horizontalLayout_7.addLayout(self.verticalLayout_5)
        self.verticalLayout_6 = QtWidgets.QVBoxLayout()
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.textEdit_9 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_9.setObjectName("textEdit_9")
        self.textEdit_9.setMinimumHeight(30)
        self.textEdit_9.setDisabled(True)
        self.verticalLayout_6.addWidget(self.textEdit_9)
        self.textEdit_10 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_10.setObjectName("textEdit_10")
        self.textEdit_10.setMinimumHeight(30)
        self.textEdit_10.setDisabled(True)
        self.verticalLayout_6.addWidget(self.textEdit_10)
        self.horizontalLayout_7.addLayout(self.verticalLayout_6)
        self.verticalLayout_14.addLayout(self.horizontalLayout_7)
        self.horizontalLayout_9 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_9.setObjectName("horizontalLayout_9")
        self.label_17 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_17.setObjectName("label_17")
        self.label_17.setMinimumWidth(150)
        self.horizontalLayout_9.addWidget(self.label_17)
        self.verticalLayout_9 = QtWidgets.QVBoxLayout()
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.label_18 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_18.setObjectName("label_18")
        self.label_18.setMinimumWidth(150)
        self.verticalLayout_9.addWidget(self.label_18)
        self.label_19 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_19.setObjectName("label_19")
        self.label_19.setMinimumWidth(150)
        self.verticalLayout_9.addWidget(self.label_19)
        self.horizontalLayout_9.addLayout(self.verticalLayout_9)
        self.verticalLayout_10 = QtWidgets.QVBoxLayout()
        self.verticalLayout_10.setObjectName("verticalLayout_10")
        self.textEdit_13 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_13.setObjectName("textEdit_13")
        self.textEdit_13.setMinimumHeight(30)
        self.textEdit_13.setDisabled(True)
        self.verticalLayout_10.addWidget(self.textEdit_13)
        self.textEdit_14 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_14.setObjectName("textEdit_14")
        self.textEdit_14.setMinimumHeight(30)
        self.textEdit_14.setDisabled(True)
        self.verticalLayout_10.addWidget(self.textEdit_14)
        self.horizontalLayout_9.addLayout(self.verticalLayout_10)
        self.verticalLayout_14.addLayout(self.horizontalLayout_9)
        self.horizontalLayout_10 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_10.setObjectName("horizontalLayout_10")
        self.label_20 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_20.setObjectName("label_20")
        self.label_20.setMinimumWidth(150)
        self.horizontalLayout_10.addWidget(self.label_20)
        self.verticalLayout_11 = QtWidgets.QVBoxLayout()
        self.verticalLayout_11.setObjectName("verticalLayout_11")
        self.label_21 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_21.setObjectName("label_21")
        self.label_21.setMinimumWidth(150)
        self.verticalLayout_11.addWidget(self.label_21)
        self.label_22 = QtWidgets.QLabel(parent=self.verticalLayoutWidget_4)
        self.label_22.setObjectName("label_22")
        self.label_22.setMinimumWidth(150)
        self.verticalLayout_11.addWidget(self.label_22)
        self.horizontalLayout_10.addLayout(self.verticalLayout_11)
        self.verticalLayout_12 = QtWidgets.QVBoxLayout()
        self.verticalLayout_12.setObjectName("verticalLayout_12")
        self.textEdit_15 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_15.setObjectName("textEdit_15")
        self.textEdit_15.setMinimumHeight(30)
        self.textEdit_15.setDisabled(True)
        self.verticalLayout_12.addWidget(self.textEdit_15)
        self.textEdit_16 = QtWidgets.QTextEdit(parent=self.verticalLayoutWidget_4)
        self.textEdit_16.setObjectName("textEdit_16")
        self.textEdit_16.setMinimumHeight(30)
        self.textEdit_16.setDisabled(True)
        self.verticalLayout_12.addWidget(self.textEdit_16)
        self.horizontalLayout_10.addLayout(self.verticalLayout_12)
        self.verticalLayout_14.addLayout(self.horizontalLayout_10)
        self.pushButton = QtWidgets.QPushButton(parent=self.verticalLayoutWidget_4)
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setMinimumHeight(50)
        self.pushButton.clicked.connect(self.generate)
        self.pushButton.clicked.connect(MainWindow.close)
        self.verticalLayout_14.addWidget(self.pushButton)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(parent=MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 524, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(parent=MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label_23.setText(_translate("MainWindow", "Report"))
        self.label_6.setText(_translate("MainWindow", "Experiment"))
        self.label.setText(_translate("MainWindow", "Reference Genome"))
        self.label_2.setText(_translate("MainWindow", "gtf Annotation"))
        self.label_8.setText(_translate("MainWindow", "Trimming"))
        self.label_9.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_10.setText(_translate("MainWindow", "Parameter/Values"))
        self.label_24.setText(_translate("MainWindow", "Data Processing"))
        self.label_25.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_26.setText(_translate("MainWindow", "Parameter/ Values"))
        self.label_27.setText(_translate("MainWindow", "Alignment"))
        self.label_28.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_29.setText(_translate("MainWindow", "Parameter/ Values"))
        self.label_17.setText(_translate("MainWindow", "Isoform Detection"))
        self.label_12.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_13.setText(_translate("MainWindow", "Parameter/ Values"))
        self.label_11.setText(_translate("MainWindow", "Post-Processing"))
        self.label_18.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_19.setText(_translate("MainWindow", "Parameter/ Values"))
        self.label_20.setText(_translate("MainWindow", "Statistics"))
        self.label_21.setText(_translate("MainWindow", "Method/ Tool"))
        self.label_22.setText(_translate("MainWindow", "Parameter/ Values"))
        self.textEdit_7.setText(_translate("MainWindow", "PoreChop"))
        self.textEdit_8.setText(_translate("MainWindow", "none"))
        self.textEdit_17.setText(_translate("MainWindow", "Seqkit"))
        self.textEdit_18.setText(_translate("MainWindow", "-m 100 (Minimum length)"))
        self.textEdit_19.setText(_translate("MainWindow", "DeSALT"))
        self.textEdit_20.setText(_translate("MainWindow", "-s 2 (Seed step), -l 14 (Seeding l-mer), -x ont1d (Read type)"))
        self.textEdit_13.setText(_translate("MainWindow", "FLAIR"))
        self.textEdit_14.setText(_translate("MainWindow", "FLAIR-correct, FLAIR-quantify -quality 4"))
        self.textEdit_9.setText(_translate("MainWindow", "Samtools"))
        self.textEdit_10.setText(_translate("MainWindow", "none"))
        self.textEdit_15.setText(_translate("MainWindow", "ZO-transformed, Bonferroni-corrected Beta Regression"))
        self.textEdit_16.setText(_translate("MainWindow", "Significance niveau: 0.05"))
        self.pushButton.setText(_translate("MainWindow", "Generate"))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    Dialog = QtWidgets.QMainWindow()
    ui = Ui_Report()
    ui.setupUi(Dialog)
    Dialog.show()
    sys.exit(app.exec())