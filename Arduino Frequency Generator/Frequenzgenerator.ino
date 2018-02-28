//--------------------------------------------------------------------------------
//Display initialisieren und Header Dateien für dessen Verwendung Laden

#include <Wire.h>                 //für Display 
#include <LiquidCrystal_I2C.h>    //für Display
LiquidCrystal_I2C lcd(0x27,20,4); //Initialisieren des Displays
//--------------------------------------------------------------------------------
//Definition der Signale

const int maxWaveform=4;      //Anzahl der Signalarten
const int maxSamplesNum=120;  //Schritte

const int waveformsTable[maxWaveform][maxSamplesNum] = {
  // Sinuswelle
  {
    0x7ff, 0x86a, 0x8d5, 0x93f, 0x9a9, 0xa11, 0xa78, 0xadd, 0xb40, 0xba1,
    0xbff, 0xc5a, 0xcb2, 0xd08, 0xd59, 0xda7, 0xdf1, 0xe36, 0xe77, 0xeb4,
    0xeec, 0xf1f, 0xf4d, 0xf77, 0xf9a, 0xfb9, 0xfd2, 0xfe5, 0xff3, 0xffc,
    0xfff, 0xffc, 0xff3, 0xfe5, 0xfd2, 0xfb9, 0xf9a, 0xf77, 0xf4d, 0xf1f,
    0xeec, 0xeb4, 0xe77, 0xe36, 0xdf1, 0xda7, 0xd59, 0xd08, 0xcb2, 0xc5a,
    0xbff, 0xba1, 0xb40, 0xadd, 0xa78, 0xa11, 0x9a9, 0x93f, 0x8d5, 0x86a,
    0x7ff, 0x794, 0x729, 0x6bf, 0x655, 0x5ed, 0x586, 0x521, 0x4be, 0x45d,
    0x3ff, 0x3a4, 0x34c, 0x2f6, 0x2a5, 0x257, 0x20d, 0x1c8, 0x187, 0x14a,
    0x112, 0xdf, 0xb1, 0x87, 0x64, 0x45, 0x2c, 0x19, 0xb, 0x2,
    0x0, 0x2, 0xb, 0x19, 0x2c, 0x45, 0x64, 0x87, 0xb1, 0xdf,
    0x112, 0x14a, 0x187, 0x1c8, 0x20d, 0x257, 0x2a5, 0x2f6, 0x34c, 0x3a4,
    0x3ff, 0x45d, 0x4be, 0x521, 0x586, 0x5ed, 0x655, 0x6bf, 0x729, 0x794
  }
  ,

  // Dreiecksignal
  {
    0x44, 0x88, 0xcc, 0x110, 0x154, 0x198, 0x1dc, 0x220, 0x264, 0x2a8,
    0x2ec, 0x330, 0x374, 0x3b8, 0x3fc, 0x440, 0x484, 0x4c8, 0x50c, 0x550,
    0x594, 0x5d8, 0x61c, 0x660, 0x6a4, 0x6e8, 0x72c, 0x770, 0x7b4, 0x7f8,
    0x83c, 0x880, 0x8c4, 0x908, 0x94c, 0x990, 0x9d4, 0xa18, 0xa5c, 0xaa0,
    0xae4, 0xb28, 0xb6c, 0xbb0, 0xbf4, 0xc38, 0xc7c, 0xcc0, 0xd04, 0xd48,
    0xd8c, 0xdd0, 0xe14, 0xe58, 0xe9c, 0xee0, 0xf24, 0xf68, 0xfac, 0xff0,
    0xfac, 0xf68, 0xf24, 0xee0, 0xe9c, 0xe58, 0xe14, 0xdd0, 0xd8c, 0xd48,
    0xd04, 0xcc0, 0xc7c, 0xc38, 0xbf4, 0xbb0, 0xb6c, 0xb28, 0xae4, 0xaa0,
    0xa5c, 0xa18, 0x9d4, 0x990, 0x94c, 0x908, 0x8c4, 0x880, 0x83c, 0x7f8,
    0x7b4, 0x770, 0x72c, 0x6e8, 0x6a4, 0x660, 0x61c, 0x5d8, 0x594, 0x550,
    0x50c, 0x4c8, 0x484, 0x440, 0x3fc, 0x3b8, 0x374, 0x330, 0x2ec, 0x2a8,
    0x264, 0x220, 0x1dc, 0x198, 0x154, 0x110, 0xcc, 0x88, 0x44, 0x0
  }
  ,

  // sägezahnsignal
  {
    0x22, 0x44, 0x66, 0x88, 0xaa, 0xcc, 0xee, 0x110, 0x132, 0x154,
    0x176, 0x198, 0x1ba, 0x1dc, 0x1fe, 0x220, 0x242, 0x264, 0x286, 0x2a8,
    0x2ca, 0x2ec, 0x30e, 0x330, 0x352, 0x374, 0x396, 0x3b8, 0x3da, 0x3fc,
    0x41e, 0x440, 0x462, 0x484, 0x4a6, 0x4c8, 0x4ea, 0x50c, 0x52e, 0x550,
    0x572, 0x594, 0x5b6, 0x5d8, 0x5fa, 0x61c, 0x63e, 0x660, 0x682, 0x6a4,
    0x6c6, 0x6e8, 0x70a, 0x72c, 0x74e, 0x770, 0x792, 0x7b4, 0x7d6, 0x7f8,
    0x81a, 0x83c, 0x85e, 0x880, 0x8a2, 0x8c4, 0x8e6, 0x908, 0x92a, 0x94c,
    0x96e, 0x990, 0x9b2, 0x9d4, 0x9f6, 0xa18, 0xa3a, 0xa5c, 0xa7e, 0xaa0,
    0xac2, 0xae4, 0xb06, 0xb28, 0xb4a, 0xb6c, 0xb8e, 0xbb0, 0xbd2, 0xbf4,
    0xc16, 0xc38, 0xc5a, 0xc7c, 0xc9e, 0xcc0, 0xce2, 0xd04, 0xd26, 0xd48,
    0xd6a, 0xd8c, 0xdae, 0xdd0, 0xdf2, 0xe14, 0xe36, 0xe58, 0xe7a, 0xe9c,
    0xebe, 0xee0, 0xf02, 0xf24, 0xf46, 0xf68, 0xf8a, 0xfac, 0xfce, 0xff0
  }
  ,

  // Rechtecksignal
  {
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff, 0xfff,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0,
    0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0
  }
};

//für 1 - 10 Hz
#define oneHzSample 1000000/maxSamplesNum   //Zeit, die gewartet wird, um 1 Hz auszugeben
#define maxHzSample oneHzSample/10          //Zeit die gewartet wird, um 10 Hz auszugeben
                                            //Die Zahl, durch die geteilt wird bestimmt
                                            //die maximale Frequenz
//--------------------------------------------------------------------------------
//Definition sonstiger Variablen

const int buttondisplay = 2;    //Pin Nummer für Display aktualisieren
const int buttonsignal = 3;     //Pin Nummer für Signalform ändern
volatile int wave0 = 0;         //Zähler für die Signalform

int i = 0;                      //Zähler der einzelnen Spannungen, um das Signal zu bilden
float freq = 0;                 //Frequenz die am Display ausgegeben wird
int sample;                     //Timer für das Warten bis zum nächsten Punkt der Lookup Tabelle
float gain;                     //Divisor für die Amplitudenhöhe

//--------------------------------------------------------------------------------
//Wird einmal beim Starten des Arduinos ausgeführt. Dabei wird das Display initialisiert, 
//Auflösungen der Ein - und Ausgänge definiert und die Ereignisse beim Drücken eines Knopfes festgelegt. 
//Die Ereignisse sind dabei in den Funktionen am Programmende zu finden

void setup() {
  analogWriteResolution(12);  //Setzt die Auflösung auf 12 bit (4096 levels)
  analogReadResolution(12);   //Setzt die Auflösung auf 12 bit (4096 levels)

  attachInterrupt(buttondisplay, displayrefresh, RISING);  //Ereignis: Knopfdruck für die Displayaktualisierung
  attachInterrupt(buttonsignal, Select, RISING);           //Ereignis: Knopfdruck für die Signalform

  lcd.init();                       //Display initialisieren
  lcd.backlight();                  //Display Hintergrundbeleuchtung einschalten
  lcd.setCursor(0, 0);              //Cursor zum Schreiben setzen (Position, Zeile)
  lcd.print("Frequenzgenerator");   //Ausgabe von Text
  lcd.setCursor(0,1);
  lcd.print("Form:     Sinus");
  lcd.setCursor(0,2);
  lcd.print("Frequenz: ");
  lcd.setCursor(0,3);
  lcd.print("Amplitude: ");
}

//--------------------------------------------------------------------------------
//Wird als Schleife ausgeführt. Dabei werden die Spannungen der beiden Potentiometer
//für die Amplitude und Frequenz eingelesen, die Spannung für den Signalausgang ausgeben und
//bis zum nächsten Schritt gewartet

void loop() {

  sample = map(analogRead(A8), 0, 4095, maxHzSample, oneHzSample);  //Einstellung der Frequenz über die sample Zeit
  gain = map(analogRead(A9), 0, 4095, 100, 10)*0.1;                 //Einstellung der Verstärkung von Stufe 1/1 bis 1/10

  analogWrite(DAC1, waveformsTable[wave0][i]/gain);  //gibt die Spannung am Ausgang DAC1 aus

  i++;                    
  if(i == maxSamplesNum)  //Setzt i zurück, wenn das Ende des Lookuptables erreicht ist
    i = 0;

  delayMicroseconds(sample);  //Wartet, um damit die Frequenz einzustellen
}

//--------------------------------------------------------------------------------
//Wird ausgeführt, wenn der Knopf zum Wechseln der Frequenz gedrückt wird
//Damit lässt sich durch die Signalformen schalten

void Select() {                        
  wave0++;          //Variabel
  if(wave0 == 4)
    wave0 = 0;
}

//--------------------------------------------------------------------------------
//Gibt die Werte auf dem Display auf und hält diese bis zum nächsten drücken. Für diese
// Aufgabe wird Zeit benötigt. Daher wurde sie aus der Schleife ausgelagert, um die maximal
//mögliche Frequenz zu erhöhen

void displayrefresh() {
  
  lcd.setCursor(10,2);
  freq = 7900.0/sample;   //Berechnet die Frequenz aus der sample Zeit; Wert experimentell ermittelt
  lcd.print(freq);
  lcd.print(" Hz   ");
  lcd.setCursor(0,3);
  lcd.print("Amplitude:    ");
  lcd.setCursor(11,3);
  lcd.print(1.09/gain);   //Gibt die Amplitude aus. +- 1.09 V ist dabei die maximale Amplitude 
  lcd.setCursor(15,3);
  lcd.print(" V ");

  lcd.setCursor(10, 1);   //Gibt die Signalform aus
     switch (wave0){                        
      case 0:
        lcd.print("Sinus     ");
      break;  
      case 1:
        lcd.print("Dreieck   ");
      break;
      case 2:
        lcd.print("Saegezahn ");
      break;
      case 3:
        lcd.print("Rechteck  ");
      break;
  }
}

