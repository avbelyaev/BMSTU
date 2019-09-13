# ML

## Usage

- set up venv
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

## Example

```bash
python editor.py
```