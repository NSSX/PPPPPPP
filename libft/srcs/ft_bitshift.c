/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   ft_bitshift.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: igomez <igomez@student.42.fr>              +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2015/01/03 19:14:21 by igomez            #+#    #+#             */
/*   Updated: 2015/01/03 19:14:26 by igomez           ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

int			ft_bitshift_r(int n, unsigned int size)
{
	int save;
	int mask;

	if (!size)
		return (n);
	mask = ((1U << size) - 1);
	save = ~mask & n;
	n &= mask;
	if (n & 1U)
		n = (1U << (size - 1)) | (n >> 1);
	else
		n = n >> 1;
	return (n | save);
}

int			ft_bitshift_l(int n, unsigned int size)
{
	int save;
	int mask;

	if (!size)
		return (n);
	mask = ((1U << size) - 1);
	save = ~mask & n;
	n &= mask;
	if (n & (1U << (n - 1)))
		n = 1U | (n << 1);
	else
		n = n << 1;
	return (n | save);
}

int			ft_bitnshift(int n, int n_bit, unsigned int size)
{
	unsigned int j;
	unsigned int k;

	if (!size)
		return (n);
	k = (unsigned int)((n_bit < 0) ? -n_bit : n_bit);
	j = 0;
	while (j++ < k)
		n = (n_bit < 0) ? ft_bitshift_l(n, size) : ft_bitshift_r(n, size);
	return (n);
}
