import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, CheckButtons


def newton_method(z, c, max_iter=100, tol=1e-6):
    """Применяет метод Ньютона для нахождения корня уравнения f(z) = z^3 + c."""
    for i in range(max_iter):
        f_z = z ** 3 + c
        f_prime_z = 3 * z ** 2

        if f_prime_z == 0:
            break

        z_next = z - f_z / f_prime_z

        if abs(z_next - z) < tol:
            return z_next, i

        z = z_next

    return z, max_iter


def compute_basins(c, grid_size=500, R=2, center=(0, 0), max_iter=100):
    """Вычисляет бассейны Ньютона в комплексной плоскости."""
    x_center, y_center = center
    x = np.linspace(x_center - R, x_center + R, grid_size)
    y = np.linspace(y_center - R, y_center + R, grid_size)
    basins = np.zeros((grid_size, grid_size))
    iterations = np.zeros((grid_size, grid_size))

    roots = []
    root_tolerance = 1e-5

    for i, real in enumerate(x):
        for j, imag in enumerate(y):
            z = complex(real, imag)
            root, iter_count = newton_method(z, c, max_iter)

            for index, r in enumerate(roots):
                if abs(root - r) < root_tolerance:
                    basins[i, j] = index
                    iterations[i, j] = iter_count
                    break
            else:
                roots.append(root)
                basins[i, j] = len(roots) - 1
                iterations[i, j] = iter_count

    return basins, iterations


def plot_basins(ax, basins, R, center=(0, 0), show_grid=True):
    """Отображает картинку с бассейнами Ньютона."""
    x_center, y_center = center
    extent = [x_center - R, x_center + R, y_center - R, y_center + R]
    ax.imshow(basins.T, extent=extent, origin='lower', cmap='nipy_spectral')
    ax.set_title("Newton's Fractals for f(z) = z^3 + c")
    ax.set_xlabel('Re(z)')
    ax.set_ylabel('Im(z)')
    if show_grid:
        ax.grid(color='white', linestyle='--', linewidth=0.5)


def update(val):
    """Обновляет изображение при изменении параметров."""
    c_real = c_slider.val
    c_imag = c_imag_slider.val
    R = R_slider.val
    show_grid = grid_check.get_status()[0]

    c = complex(c_real, c_imag)
    basins, _ = compute_basins(c, grid_size=500, R=R, center=center, max_iter=100)

    ax.clear()  # Очищаем старое изображение
    plot_basins(ax, basins, R, center=center, show_grid=show_grid)
    fig.canvas.draw_idle()


def on_scroll(event):
    """Масштабирует изображение при прокрутке."""
    global R
    if event.button == 'up':
        R /= 1.2  # Уменьшаем область для увеличения масштаба
    elif event.button == 'down':
        R *= 1.2  # Увеличиваем область для уменьшения масштаба
    R_slider.set_val(R)


def on_drag(event):
    """Перемещает изображение при перетаскивании мышью."""
    global center
    if event.button == 1 and event.inaxes:  # Левая кнопка мыши
        dx = (event.xdata - center[0]) / 20
        dy = (event.ydata - center[1]) / 20
        center = (center[0] - dx, center[1] - dy)
        update(None)


def main():
    global c_slider, c_imag_slider, R_slider, grid_check, fig, ax, R, center

    # Начальные параметры
    c = complex(1, 1)
    R = 2
    center = (0, 0)
    grid_size = 500

    # Создание начального графика
    fig, ax = plt.subplots()
    plt.subplots_adjust(left=0.25, bottom=0.35)
    basins, _ = compute_basins(c, grid_size, R, center)
    plot_basins(ax, basins, R, center)

    # Добавление слайдеров для изменения параметров
    axcolor = 'lightgoldenrodyellow'
    ax_c_real = plt.axes([0.25, 0.2, 0.65, 0.03], facecolor=axcolor)
    ax_c_imag = plt.axes([0.25, 0.15, 0.65, 0.03], facecolor=axcolor)
    ax_R = plt.axes([0.25, 0.1, 0.65, 0.03], facecolor=axcolor)

    c_slider = Slider(ax_c_real, 'Re(c)', -2.0, 2.0, valinit=c.real)
    c_imag_slider = Slider(ax_c_imag, 'Im(c)', -2.0, 2.0, valinit=c.imag)
    R_slider = Slider(ax_R, 'R', 0.1, 4.0, valinit=R)

    # Чекбокс для отображения сетки
    ax_checkbox = plt.axes([0.025, 0.5, 0.15, 0.15], facecolor=axcolor)
    grid_check = CheckButtons(ax_checkbox, ['Grid'], [True])

    # Кнопка обновления
    ax_button = plt.axes([0.8, 0.025, 0.1, 0.04])
    button = Button(ax_button, 'Update', color=axcolor, hovercolor='0.975')

    # Привязка обновлений
    c_slider.on_changed(update)
    c_imag_slider.on_changed(update)
    R_slider.on_changed(update)
    grid_check.on_clicked(update)
    button.on_clicked(update)

    # Привязка событий прокрутки и перетаскивания
    fig.canvas.mpl_connect('scroll_event', on_scroll)
    fig.canvas.mpl_connect('motion_notify_event', on_drag)

    plt.show()


if __name__ == "__main__":
    main()
