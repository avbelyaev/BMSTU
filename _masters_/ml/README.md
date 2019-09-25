# ML

## Image editor. lab 1
#### Usage

- set up & activate venv
```bash
source venv/bin/activate

pip install -r requirements.txt
```

- use `ImgEditor` class in builder-like way:

```python
ImgEditor('source.png') \
    .rotate(RotationPoint.CENTER, angle=10) \
    .shrink(x_factor=0.9) \
    .translate(y_offset=20) \
    .save('dest.png')
```

#### Example

```bash
cd image-edit
python editor.py
```


## Square detector. lab 2

Для выбранных изображений (см почту) реализовать распознавание границ геометрических фигур. 
Распознанные границы нанести на рисунок (другим цветом).

#### Usage

- set up & activate venv
```bash
source venv/bin/activate

pip install -r requirements.txt
```

- run lab
```bash
cd square-contours
python contours.py
```
