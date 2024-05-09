# Модель взаимодействия солнечного ветра с атмосферой планеты
## Выполнена с использованием библиотеки **FEniCS**

Микропроекты находятся в директориях `atmosphere/` и `magnetic/`. Каждый из них решают подзадачу, необходимую для реализации полноценной модели (макропроекта).

Макропроект расположен в директории `macro/`.

Сборка проекта в любой из директорий осуществляется стандартным образом:
```
mkdir build
cd build
cmake ..
make
```

После запуска результаты вычислений будут сохраняться в директорию `results/`, ***которая очищается при повторном запуске***.

Также, файл `macro/Constants.hpp` содержит параметры модели макропроекта, с помощью которых ее можно легко настроить под конкретную задачу.

### Авторы
*Забавное* распределение давления в атмосфере: [костя ломанов](t.me/gazasd) <br> Разваливающееся магнитное поле: [киша малашников](t.me/laythe)