# form-SANC
немного исходников от SANC.jinr.ru

# внесенные изменеия
* переписал Declar.h в Declare.h, Parameters.h, Misc.h, добавил туда #ifndef `$FileH' #$FileH = 1 ;
* Перегруппировал a2b.prc
* Оставил только необходимое в FeynmanRules.prc
* Удалил QCDAlgeba.prc, и закомментировал ее вызовы
* Сделал возможность быстро переключаться между соглашениями SANC и Пескина-Наумова (#define PeskinNaumov "1")
 *	в пропагаторах
 *	в уравнении Дирака
 *	в Hermitian-е
 *	в спиновых суммах
* Добавил преобразование p^2 -> m^2 в Convert.prc (#$mom2mas = 1;)

#навигация по файлам

* Declare.h - термины, в которых изначально записывается матричный элемент
* Parameters.h - массы, заряды и другие параметры частиц
* Globals.prc - массы, заряды и другие параметры частиц, записанные в массивы (таблицы на языке форма)
* Misc.h - переменные и функции, которые по разному используются в разных местах

---
* a2b.prc - преобразование одних констант в другие по фиксированным формулам
* Convert.prc - преобразование нумерованных величин (масс например) частиц в именованные
* GammaLeft.prc/GammaRight.prc - протаскивание gamma5, gamma6 и gamma7 влево/вправо 
* SetFlags.prc - похоже на какую-то заглушку, которая постоянно используется
* Stop.prc - останов в случае ошибки

---
* FeynmanRules.prc - раскрываем термины, в которых изначально записывается матричный элемент
* DiracEquation.prc - некоторое упрощение матричного элемента при помощи ур-я Дирака
* Hermitian.prc - вычисляем комплексное сопряжение
* MakeAmpSquare.prc - заменяем соотв. индексы, и умножаем исдный на результат Hermitian-а 
* Trace.prc - раскрываем спиновые суммы и считаем trace

---
* ComptonBorn.frm - здесь уже кинематика
* BhabhaBorn.frm - здесь уже кинематика

#todo
* выяснить, почему исходники Рената для комптона выдают неправильный результат
* выяснить, почему модифицированные исходники в режиме SANC для комптона выдают неправильный результат
