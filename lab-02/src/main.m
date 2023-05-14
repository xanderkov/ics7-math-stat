function main()
	pkg load statistics

	X = dlmread("input.txt", ",");

    gamma = 0.8;
    n = length(X);
    % Точечная оценка мат. ожидания
    mu = mean(X);
    % Точечная оценка дисперсии
    s2 = var(X);

    % Нижняя граница доверительного интервала для мат. ожидания
    mu_low = get_mu_low(n, mu, s2, gamma);
    % Верхняя граница доверительного интервала для мат. ожидания
    mu_high = get_mu_high(n, mu, s2, gamma);
    % Нижняя граница доверительного интервала для дисперсии
    s2_low = get_s2_low(n, s2, gamma);
    % Верхняя граница доверительного интервала для дисперсии
    s2_high = get_s2_high(n, s2, gamma);

    % Вывод полученных ранее значений
    fprintf('a) Точечная оценка математического ожидания = %.3f\n', mu);
    fprintf('   Точечная оценка дисперсии = %.3f\n', s2);
    fprintf('б) Нижняя граница доверительного интервала для математического ожидания = %.3f\n', mu_low);
    fprintf('   Верхняя граница доверительного интервала для математического ожидания = %.3f\n', mu_high);
    fprintf('в) Нижняя граница доверительного интервала для дисперсии = %.3f\n', s2_low);
    fprintf('   Верхняя граница доверительного интервала для дисперсии = %.3f\n', s2_high);

    % Создание массивов точечных оценок

    n_start = 1;

    mus = zeros(n_start, n);
    s2s = zeros(n_start, n);
    % Создание массивов границ доверительных интервалов
    mus_low = zeros(n_start, n);
    mus_high = zeros(n_start, n);
    s2s_low = zeros(n_start, n);
    s2s_high = zeros(n_start, n);

    for i = n_start:n
        mu = mean(X(n_start:i));
        s2 = var(X(n_start:i));
        % Точечная оценка матожидания
        mus(i) = mu;
        % Точечная оценка дисперсии
        s2s(i) = s2;
        % Нижняя граница доверительного интервала для матожидания
        mus_low(i) = get_mu_low(i, mu, s2, gamma);
        % Верхняя граница доверительного интервала для матожидания
        mus_high(i) = get_mu_high(i, mu, s2, gamma);
        % Нижняя граница доверительного интервала для дисперсии
        s2s_low(i) = get_s2_low(i, s2, gamma);
        % Верхняя граница доверительного интервала для дисперсии
        s2s_high(i) = get_s2_high(i, s2, gamma);
    end

    % Построение графиков
    plot(n_start:n, [(zeros(n_start, n) + mu)', mus', mus_low', mus_high']);
    grid on;
    xlabel('n');
    ylabel('y');

    legend('точечная оценка \mu(x_N)', 'точечная оценка \mu(x_n)', 'верхняя граница ДО \mu', ...
        'нижняя граница ДО \mu', 'Interpreter', 'tex');
    set(gca, "ylim", [-5, -4]);
    set(gca, "xlim", [10, 150]);
    print -djpg g-1.jpg

    figure;
    plot(n_start:n, [(zeros(1, n) + s2)', s2s', s2s_low', s2s_high']);
    grid on;
    xlabel('n');
    ylabel('z');
    legend('точечная оценка S^2(x_N)', 'точечная оценка S^2(x_n)', 'верхняя граница ДО \sigma', ...
            'нижняя граница ДО \sigma', 'Interpreter', 'tex');
    set(gca, "ylim",[0, 2]);
    set(gca, "xlim", [10, 150]);
    print -djpg g-2.jpg
end

function mu_low = get_mu_low(n, mu, s2, gamma)
    mu_low = mu - sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function mu_high = get_mu_high(n, mu, s2, gamma)
    mu_high = mu + sqrt(s2) * tinv((1 + gamma) / 2, n - 1) / sqrt(n);
end

function s2_low = get_s2_low(n, s2, gamma)
    s2_low = ((n - 1) * s2) / chi2inv((1 + gamma) / 2, n - 1);
end

function s2_high = get_s2_high(n, s2, gamma)
    s2_high = ((n - 1) * s2) / chi2inv((1 - gamma) / 2, n - 1);
end
