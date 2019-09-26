# ML

Prerequirements

- set up & activate venv
```bash
source venv/bin/activate
pip install -r requirements.txt
```

## Image editor. lab 1

Аугументация.

Реализовать программу (python или C++), на вход которой подаётся изображения, а на выходе должно быть несколько изображений:

- Повороты этого изображения на некоторый угол
- Смещение этого изображения по оси x и по оси y
- Сжатие и растяжение
- Параллельный перенос
- Симметрия

#### Usage

```python
ImgEditor('source.png') \
    .rotate(RotationPoint.CENTER, angle=10) \
    .shrink(x_factor=0.9) \
    .translate(y_offset=20) \
    .save('dest.png')
```

#### Example

```bash
cd image_edit
python editor.py
```


## Square detector. lab 2

Для выбранных изображений (см почту) реализовать распознавание границ геометрических фигур. 
Распознанные границы нанести на рисунок (другим цветом).

#### Usage

```bash
cd square_contours
python contours.py
```


## Blur detector. lab 3

Реализовать алгоритмы распознавания размытого изображения и региона

- Tenengrad (TENG)
- Normalized Gray Level Variance (GLVN)

возвращаемое значение – метка изображение размыто или нет. 
Для региона указать координаты размытого региона на изображении.
