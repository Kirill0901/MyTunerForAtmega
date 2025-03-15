#ifndef LCD_H_
#define LCD_H_

#define LCDDATAPORT        PORTD  // ���� � ����,
#define LCDDATADDR         DDRD   // � ������� ����������
#define LCDDATAPIN         PIND   // ������� D4-D7.
#define LCD_D4             0
#define LCD_D5             1
#define LCD_D6             2
#define LCD_D7             3

#define LCDCONTROLPORT     PORTD  // ���� � ����,
#define LCDCONTROLDDR      DDRD   // � ������� ����������
#define LCD_RS             6      // ������� RS, RW � E.
#define LCD_RW             5
#define LCD_E              4

#define LCD_STROBEDELAY_US 5      // �������� ������

#define LCD_COMMAND        0
#define LCD_DATA           1

#define LCD_CURSOR_OFF     0
#define LCD_CURSOR_ON      2
#define LCD_CURSOR_BLINK   3

#define LCD_DISPLAY_OFF    0
#define LCD_DISPLAY_ON     4

#define LCD_SCROLL_LEFT    0
#define LCD_SCROLL_RIGHT   4

#define LCD_STROBDOWN      0
#define LCD_STROBUP        1

#define DELAY              1
 
void lcdSendNibble(char byte, char state);
char lcdGetNibble(char state);
char lcdRawGetByte(char state);
void lcdRawSendByte(char byte, char state);
char lcdIsBusy(void);
void lcdInit(void);
void lcdSetCursor(char cursor);
void lcdSetDisplay(char state);
void lcdClear(void);
void lcdGotoXY(char str, char col);
void lcdDisplayScroll(char pos, char dir);
void lcdPuts(char *str);
void lcdPutsf(char *str);
void lcdPutse(uint8_t *str);
void lcdLoadCharacter(char code, char *pattern);
void lcdLoadCharacterf(char code, char *pattern);
void lcdLoadCharactere(char code, char *pattern);

#endif /* LCD_H_ */